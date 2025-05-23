#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <string>
#include <map>

#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_parallel.h> // parallel N_Vector types, fcts., and macros

#ifdef SUNDIALS2
#include <cvode/cvode_spgmr.h>
#elif SUNDIALS3
#include <cvode/cvode_spils.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#elif SUNDIALS4
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <cvode/cvode_bbdpre.h>
#endif

#include "utilities.h"
#include <utilities/file_utilities.h>

#include <mpi.h>

#include "flame_params.h"
#include "cvode_functions.h"
#include "set_initial_conditions.h"

using zerork::getHighResolutionTime;

const bool RUN_DEBUG=false;
const int NUM_STDOUT_PARAMS = 19; // number of non species parameters to write
                                  // to standard out

static int check_flag(void *flagvalue, const char *funcname, int opt);

static double FindMaximumParallel(const size_t num_points,
				  const double x[],
				  const size_t x_stride,
				  const double f[],
				  const size_t f_stride,
				  const bool use_quadratic,
				  double *x_at_max);

static void WriteFieldParallel(double t,
			       const double state[],
			       const FlameParams &params);

static double GetMixtureMolecularMass(const int grid_id,
                                      const double t,
                                      const double state[],
                                      const FlameParams &params);
static double minSumMassFractions(const double state[],
                                  const FlameParams &params);
static double maxSumMassFractions(const double state[],
                                  const FlameParams &params);
static double minVelocity(const double state[],
			  const FlameParams &params);
static double maxVelocity(const double state[],
			  const FlameParams &params);
static double GetDensityFromYTP(const int grid_id,
                                const double t,
                                const double state[],
                                const FlameParams &params);

static int GetFuelSpeciesId(const FlameParams &params,
                            std::vector<int> *species_id);

static double InletFuelFraction(const FlameParams &params);


static int GetTrackMaxStateId(const FlameParams &params,
                              std::vector<int> *state_id);
static int GetStateMaxima(const std::vector<int> &state_id,
                          const double state[],
                          const FlameParams &params,
                          std::vector<double> *max_value,
                          std::vector<double> *max_position);

static double StoichiometricMixtureFraction(const FlameParams &params);

static int SootOutput(const FlameParams &params,  const double state[]);

static double Min(double a, double b)
{
  return ((a < b) ? a : b);
}

int main(int argc, char *argv[])
{
  double clock_time = getHighResolutionTime();
  double setup_time, loop_time;

  if(argc < 2) {
    printf("# ERROR: Incorrect command line usage.\n");
    printf("#        use %s <input parameters>\n",argv[0]);
    exit(-1);
  }

  // MPI
  MPI_Comm comm;
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;

  // CVode memory pointer and linear solver
  void *cvode_ptr = NULL;
#ifdef SUNDIALS3
  SUNLinearSolver LS;
  LS = NULL;
#elif SUNDIALS4
  SUNLinearSolver LS;
  SUNNonlinearSolver NLS;
  LS = NULL;
  NLS = NULL;
#endif

  // Initialize flame params
  FlameParams flame_params(argv[1],comm);

  // Initialize state variables vector and pointer
  N_Vector flame_state;
  flame_state = NULL;
  double *flame_state_ptr;
  flame_state_ptr = NULL;

  // Get constants from flame params
  const int num_grid_points = flame_params.z_.size();
  const int num_local_points = flame_params.num_local_points_;
  const int num_states = flame_params.reactor_->GetNumStates();
  const int num_species = flame_params.reactor_->GetNumSpecies();
  long int num_local_states = num_local_points*num_states;
  long int total_states = num_grid_points*num_states;
  double max_time = flame_params.parser_->max_time();
  const double step_dt = flame_params.parser_->stats_dt();
  const double ref_temperature = flame_params.ref_temperature_;
  const double dz = flame_params.dz_[(int)num_grid_points/2]; //should be min(dz_)
  int my_pe = flame_params.my_pe_;
  int npes = flame_params.npes_;

  // Declare variables
  int flag=0;
  int num_prints = 0;
  double current_time = 0.0;
  double time_offset;
  double next_time;
  double max_temperature_jump, z_max_temperature_jump;
  double fuel_density, fuel_molecular_mass;
  double oxidizer_density, oxidizer_molecular_mass;
  double stoichiometric_mixture_fraction;
  double min_sum_mass_fraction, max_sum_mass_fraction;
  double min_velocity, max_velocity;

  // CVode integrator Stats
  long int nsteps, nfevals, nlinsetups, netfails;
  int qlast, qcur;
  double hinused, hlast, hcur, tcur;
  std::vector<double> temperature_jump;
  std::vector<int> track_max_state_id;
  std::vector<double> state_maxima, state_maxima_positions;


  // allocate CVode data structures
  flame_state          = N_VNew_Parallel(comm, num_local_states, total_states);
  flame_state_ptr      = NV_DATA_P(flame_state);

  // Initialize state vector
  time_offset = 0.0;
  SetInitialComposition(flame_params, flame_state_ptr, &time_offset);
  max_time -= time_offset;

  temperature_jump.assign(num_local_points,0.0);
  GetTrackMaxStateId(flame_params, &track_max_state_id);

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration */
#ifdef SUNDIALS4
  cvode_ptr = CVodeCreate(CV_BDF);
#else
  cvode_ptr = CVodeCreate(CV_BDF, CV_NEWTON);
#endif
  if (check_flag((void *)cvode_ptr, "CVodeCreate", 0)) exit(-1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag=CVodeInit(cvode_ptr,
                 ConstPressureFlame,
                 0.0,
                 flame_state);
  if (check_flag(&flag, "CVodeInit", 1)) exit(-1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  flag = CVodeSStolerances(cvode_ptr,
                           flame_params.parser_->rel_tol(),
                           flame_params.parser_->abs_tol());
  if (check_flag(&flag, "CVodeSStolerances", 1)) exit(-1);

  /* Set the pointer to user-defined data */
  flag = CVodeSetUserData(cvode_ptr,
                          &flame_params);
  if(check_flag(&flag, "CVodeSetUserData", 1)) exit(-1);

#ifdef SUNDIALS2
  /* Setup the linear solver method */
  flag = CVSpgmr(cvode_ptr, PREC_LEFT, 5);
  if(check_flag(&flag, "CVSpgmr", 1)) exit(-1);

  flag = CVSpilsSetGSType(cvode_ptr, MODIFIED_GS);
  if(check_flag(&flag, "CVSpilsSetGSType", 1)) exit(-1);

   /* Set the preconditioner setup and solve functions */
  flag = CVSpilsSetPreconditioner(cvode_ptr,
                                  ReactorPreconditionerChemistrySetup,
                                  ReactorPreconditionerChemistrySolve);
  if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) exit(-1);

#elif SUNDIALS3
  /* Setup the linear solver method */
  // TODO: add SPGMR solver parameters to input file for user control
  LS = SUNSPGMR(flame_state,
                PREC_LEFT,
                5); //max krylov dimension

  flag = CVSpilsSetLinearSolver(cvode_ptr, LS);
  if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) exit(-1);

  flag = SUNSPGMRSetGSType(LS, MODIFIED_GS);
  if(check_flag(&flag, "SUNSPGMRSetGSType", 1)) exit(-1);

  /* Set the preconditioner setup and solve functions */
  flag = CVSpilsSetPreconditioner(cvode_ptr,
                                  ReactorPreconditionerChemistrySetup,
                                  ReactorPreconditionerChemistrySolve);
  if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) exit(-1);
#elif SUNDIALS4
  /* Setup the nonlinear solver method */
  NLS = SUNNonlinSol_Newton(flame_state);
  flag = CVodeSetNonlinearSolver(cvode_ptr, NLS);
  if(check_flag(&flag, "CVodeSetNonlinearSolver", 1)) exit(-1);

  /* Setup the linear solver method */
  // TODO: add SPGMR solver parameters to input file for user control
  LS = SUNLinSol_SPGMR(flame_state,
                       PREC_LEFT,
                       5); //max krylov dimension

  flag = SUNSPGMRSetGSType(LS, MODIFIED_GS);
  if(check_flag(&flag, "SUNSPGMRSetGSType", 1)) exit(-1);

  flag = CVodeSetLinearSolver(cvode_ptr, LS, NULL);
  if(check_flag(&flag, "CVodeSetLinearSolver", 1)) exit(-1);

  // Set the preconditioner setup and solve functions
  flag = CVodeSetPreconditioner(cvode_ptr,
                                ReactorPreconditionerChemistrySetup,
                                ReactorPreconditionerChemistrySolve);
  if(check_flag(&flag, "CVodeSetPreconditioner", 1)) exit(-1);

#endif

  /* Set the maximum number of internal steps per CVode call and the maximum
   * allowable internal steps. */
  flag = CVodeSetMaxNumSteps(cvode_ptr,
                             flame_params.parser_->max_num_steps());
  if(check_flag(&flag, "CVodeSetMaxNumSteps", 1)) exit(-1);

  flag = CVodeSetMaxStep(cvode_ptr,
                         flame_params.parser_->max_internal_dt());
  if(check_flag(&flag, "CVodeSetMaxStep", 1)) exit(-1);

  // Write initial data file
  if (time_offset == 0.0) {
    WriteFieldParallel(time_offset,
		       &flame_state_ptr[0],
		       flame_params);
    SootOutput(flame_params,&flame_state_ptr[0]);
  }

  // Report simulation info
  if(my_pe == 0) {
    fuel_molecular_mass = GetMixtureMolecularMass(0, 0.0,
	   &flame_params.fuel_mass_fractions_[0],
	   flame_params);
    oxidizer_molecular_mass = GetMixtureMolecularMass(0, 0.0,
                                                  &flame_params.oxidizer_mass_fractions_[0],
                                                  flame_params);
    fuel_density = flame_params.parser_->pressure()*fuel_molecular_mass/
      (flame_params.fuel_temperature_*flame_params.ref_temperature_*
       flame_params.reactor_->GetGasConstant());
    oxidizer_density = flame_params.parser_->pressure()*oxidizer_molecular_mass/
      (flame_params.oxidizer_temperature_*flame_params.ref_temperature_*
       flame_params.reactor_->GetGasConstant());
    stoichiometric_mixture_fraction = StoichiometricMixtureFraction(flame_params);
    printf("# Number of states     : %d\n",num_states);
    printf("# Number of grid points: %d\n", num_grid_points);
    printf("# Stoichiometric Mixture fraction: %.18g\n",
           stoichiometric_mixture_fraction);
    printf("# Fuel BC temperature  [K]: %.18g\n",
           flame_params.fuel_temperature_*flame_params.ref_temperature_);
    printf("# Fuel BC density [kg/m^3]: %.18g\n", fuel_density);
    printf("# Oxidizer BC density [kg/m^3]: %.18g\n", oxidizer_density);
    printf("# Initial upstream temperature  (next to inlet)     [K]: %.18g\n",
           flame_state_ptr[num_species+1]*ref_temperature);
    printf("# Initial upstream density      (next to inlet)[kg/m^3]: %.18g\n",
           1.0/flame_state_ptr[num_species]);
    printf("# Initial upstream avg velocity (next to inlet)   [m/s]: %.18g\n",
           flame_params.mass_flux_ext_[2]*flame_state_ptr[num_species]);
    printf("#------------------------------------------------------------------------------\n");
    printf("# Column  1: [#] time steps printed\n");
    printf("# Column  2: [s] ODE system time\n");
    printf("# Column  3: [m/s] fuel burning rate/flame speed\n");
    printf("# Column  4: [m] flame thickness defined as T_burnt-T_unburnt\n"
	   "#            divided by maximum temperature gradient\n");
    printf("# Column  5: [K] maximum jump (T-T_wall) in the domain\n");
    printf("# Column  6: [m] location of the stagnation plane\n");
    printf("# Column  7: [#] number of steps taken by CVode\n");
    printf("# Column  8: [#] number of calls to the user's time derivative (RHS) function\n");
    printf("# Column  9: [#] number of calls to the linear solver setup function\n");
    printf("# Column 10: [#] number of error test failures\n");
    printf("# Column 11: [#] current method order of time integrator\n");
    printf("# Column 12: [s] current size of internal time step\n");
    printf("# Column 13: [s] characteristic time step for Courant number control dx/max(|velocity|) \n");
    printf("# Column 14: [s] characteristic time step for diffusion stability control 0.5*dx*dx/max(thermal diffusivity)\n");
    //printf("# Column 15: Mass flux at last grid point\n");
    printf("# Column 15: [1/s] strain rate \n");
    printf("# Column 16: min sum of mass fractions\n");
    printf("# Column 17: max sum of mass fractions\n");
    printf("# Column 18: [m/s] min velocity\n");
    printf("# Column 19: [m/s] max velocity\n");
    for(int j=0; j<(int)track_max_state_id.size(); ++j) {
      printf("# Column %2d: maximum %s in the domain\n",
	     NUM_STDOUT_PARAMS+1+2*j,
	     flame_params.reactor_->GetNameOfStateId(track_max_state_id[j]));
      printf("# Column %2d: [m] location of the maximum %s\n",
	     NUM_STDOUT_PARAMS+2+2*j,
	     flame_params.reactor_->GetNameOfStateId(track_max_state_id[j]));
    }
  }// if(my_pe==0)

  // integrate ODE system
  next_time = Min(step_dt, max_time);

  // Save time for CPU time evaluation
  setup_time = getHighResolutionTime() - clock_time;
  clock_time = getHighResolutionTime();

  // While simulation loop
  while(current_time < max_time) {
    flag = CVode(cvode_ptr,
                 next_time,
                 flame_state,
                 &current_time,
                 CV_NORMAL);

    if(flag == CV_SUCCESS || current_time >= next_time) {
      next_time += step_dt;
      ++num_prints;

      // Compute T
      for(int j=0; j<num_local_points; ++j) {
        temperature_jump[j] =
          flame_state_ptr[j*num_states+num_species+1]*flame_params.ref_temperature_;
      }
      // Find maximum value over domain
      max_temperature_jump = FindMaximumParallel(num_local_points,
						 &flame_params.z_[0],
						 1,
						 &temperature_jump[0],
						 1,
						 true, // use quadratic
						 &z_max_temperature_jump);

      // Get min/max sum(Yi) and velocity
      min_sum_mass_fraction = minSumMassFractions(flame_state_ptr,flame_params);
      max_sum_mass_fraction = maxSumMassFractions(flame_state_ptr,flame_params);
      min_velocity = minVelocity(flame_state_ptr,flame_params);
      max_velocity = maxVelocity(flame_state_ptr,flame_params);

      // Get max and locations of tracked variables
      GetStateMaxima(track_max_state_id,
                     flame_state_ptr,
                     flame_params,
                     &state_maxima,
                     &state_maxima_positions);

      // Get integrator stats
      CVodeGetIntegratorStats(cvode_ptr, &nsteps, &nfevals,
                              &nlinsetups, &netfails, &qlast, &qcur,
                              &hinused, &hlast, &hcur, &tcur);

      double mdotlast;
      if(my_pe == npes-1)
        mdotlast = flame_params.mass_flux_[num_local_points-1];
      MPI_Bcast(&mdotlast, 1, MPI_DOUBLE, npes-1, comm);

      // Print to screen/log
      if(my_pe == 0) { //only root prints
	printf("%7d  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %6d  %6d  %6d  %6d  %6d  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
	       num_prints,
	       current_time+time_offset,
	       flame_params.flame_speed_,
	       flame_params.flame_thickness_,
	       max_temperature_jump,
	       flame_params.stagnation_plane_,//z_max_temperature_jump,
	       (int)nsteps,
	       (int)nfevals,
	       (int)nlinsetups,
	       (int)netfails,
	       qcur,
	       hcur,
	       dz/flame_params.max_velocity_,
	       0.5*dz*dz/flame_params.max_thermal_diffusivity_,
               //mdotlast,
               flame_params.strain_rate_,
	       min_sum_mass_fraction,
	       max_sum_mass_fraction,
	       min_velocity,
	       max_velocity);
	for(size_t j=0; j<track_max_state_id.size(); ++j) {
	  printf("  %14.7e  %14.7e",state_maxima[j], state_maxima_positions[j]);
	}
	printf("\n");
      } // if(my_pe==0)

      // Write data file
      if(num_prints % flame_params.parser_->field_dt_multiplier() == 0) {
	WriteFieldParallel(current_time+time_offset,
			   &flame_state_ptr[0],
			   flame_params);
        SootOutput(flame_params,&flame_state_ptr[0]);
      }

    } else {
      printf("# ERROR: In time marching loop,\n");
      printf("#        CVode returned flag = %d.\n", flag);
      printf("#        Halting integration.\n");
      break;
    }

  }
  loop_time = getHighResolutionTime() - clock_time;

  N_VDestroy_Parallel(flame_state);
  if(cvode_ptr != NULL) {
    CVodeFree(&cvode_ptr);
  }
#if defined SUNDIALS3 || defined SUNDIALS4
  SUNLinSolFree(LS);
#endif
#if defined SUNDIALS4
  SUNNonlinSolFree(NLS);
#endif
  if(my_pe == 0) {
    printf("# Simulation setup time   [s]: %12.5e\n",setup_time);
    printf("# Time in integrator loop [s]: %12.5e\n",loop_time);
  }
  flame_params.logger_->PrintF(
    "# Simulation setup time   [s]: %12.5e\n",setup_time);
  flame_params.logger_->PrintF(
    "# Time in integrator loop [s]: %12.5e\n",loop_time);

  MPI_Finalize();
  return 0;
}

// Find maximum and its location in parallel
static double FindMaximumParallel(const size_t num_points,
                          const double x[],
                          const size_t x_stride,
                          const double f[],
                          const size_t f_stride,
                          const bool use_quadratic,
                          double *x_at_max)
{
  int myrank;
  struct {
    double value;
    int index;
  } in, out;

  // Compute local maximum
  in.value = f[0];
  in.index = 0;
  for(int j=1; j<(int)num_points; ++j) {
    if(in.value < f[j*f_stride]) {
      in.value = f[j*f_stride];
      in.index = j;
    }
  }

  // Compute global maximum
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  in.index += myrank*num_points;

  MPI_Allreduce(&in,&out,1,MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

  *x_at_max = x[out.index];

  return out.value;
}

// Write field to binary file for restart/postprocessing
static void WriteFieldParallel(double t,
			       const double state[],
			       const FlameParams &params)
{
  int num_grid_points = (int)params.z_.size();
  const int num_grid_points_ext = (int)params.z_.size() + 2;
  int num_local_points = (int)params.num_local_points_;
  int num_reactor_states = params.reactor_->GetNumStates();
  int num_species = params.reactor_->GetNumSpecies();
  int num_steps = params.reactor_->GetNumSteps();
  int my_pe, npes;
  char filename[32], *basename;
  bool dump_mole_fractions = params.parser_->write_mole_fractions_to_field_files();
  bool dump_step_rates = params.parser_->write_step_rates_to_field_files();

  int disp;
  std::vector<double> buffer(num_local_points+2, 0.0);
  std::vector<char> charbuf(64,0);

  MPI_File output_file;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_pe);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  int global_size = (num_grid_points+2)*sizeof(double);
  int field_count = num_local_points;
  if(my_pe == 0) field_count +=1; //left BC
  if(my_pe == npes-1) field_count +=1; //right BC
  std::vector<int> field_index(npes,0);
  for(int j = 1; j<npes; ++j) {
    field_index[j] = field_index[j-1]+num_local_points;
    if(j==1) field_index[j] += 1; //left BC
  }
  
  basename = "data_";
  sprintf(filename, "%s%f", basename, t);

  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_RDWR,
		MPI_INFO_NULL, &output_file);

  int header_size = 0;
  // Write header
  if(my_pe == 0) {
    MPI_File_write(output_file, &num_grid_points_ext, 1, MPI_INT, MPI_STATUS_IGNORE);//num points
    int num_vars = num_reactor_states; //excludes mass flux
    if(dump_mole_fractions) num_vars += num_species;
    if(dump_step_rates) num_vars += num_steps + num_species;
    MPI_File_write(output_file, &num_vars, 1, MPI_INT, MPI_STATUS_IGNORE);//num variables
    MPI_File_write(output_file, &t, 1, MPI_DOUBLE, MPI_STATUS_IGNORE); //time
    for(int j=0; j<num_reactor_states; ++j) {
      std::string state_name = params.reactor_->GetNameOfStateId(j);
      strcpy(&charbuf[0], state_name.c_str());
      MPI_File_write(output_file, &charbuf[0], 64, MPI_CHAR, MPI_STATUS_IGNORE); //state name
    }
    header_size = 2*sizeof(int) + sizeof(double) + num_vars*sizeof(char)*64;
    if(dump_mole_fractions) {
      for(int j=0; j<num_reactor_states; ++j) {
        std::string state_name = params.reactor_->GetNameOfStateId(j);
        size_t found = state_name.find("MassFraction");
        if(found != std::string::npos) {
          state_name.replace(0, 4, "Mole");
          strcpy(&charbuf[0], state_name.c_str());
          MPI_File_write(output_file, &charbuf[0], 64, MPI_CHAR, MPI_STATUS_IGNORE); //state name
        }
      }
    }
    if(dump_step_rates) {
      for(int j=0; j<num_steps; ++j) {
	std::string step_name;
	params.reactor_->GetMechanism()->getReactionNameDirOfStep(j, &step_name);
	snprintf(&charbuf[0], 64, "%s ROP", step_name.c_str());
        MPI_File_write(output_file, &charbuf[0], 64, MPI_CHAR, MPI_STATUS_IGNORE); //state name
      }
      for(int j=0; j<num_species; ++j) {
        std::string state_name = params.reactor_->GetNameOfStateId(j);
        size_t found = state_name.find("MassFraction");
        if(found != std::string::npos) {
          state_name.replace(0, 12, "NetProductionRate");
	  snprintf(&charbuf[0], 64, "%s", state_name.c_str());
          MPI_File_write(output_file, &charbuf[0], 64, MPI_CHAR, MPI_STATUS_IGNORE); //state name
	}
      }
    }
    MPI_File_seek(output_file, 0, MPI_SEEK_SET);
  }
  MPI_Bcast(&header_size,1,MPI_INT,0,MPI_COMM_WORLD);

  // Write data for each variable
  for(int j=0; j<num_reactor_states; ++j) {
    int offset = 0;
    if(my_pe == 0) {
      offset+=1;
      // Write left (fuel) BC data
      if(j==num_species){
        if(params.flame_type_ == 0) {
          buffer[0] = params.fuel_relative_volume_;
        } else {
          buffer[0] = params.inlet_relative_volume_;
        }
      } else if(j==num_species+1){
        buffer[0] = params.fuel_temperature_*params.ref_temperature_;
      } else if (j==num_species+2) {
        buffer[0] = 0.0;
      } else if (j==num_species+3) {
        buffer[0] = params.P_left_*params.ref_momentum_;
      } else {
        if(params.flame_type_ == 0) {
          buffer[0] = params.fuel_mass_fractions_[j];
        } else {
          buffer[0] = params.inlet_mass_fractions_[j];
        }
      }
    }

    // Write interior data
    for (int k=0; k<num_local_points; ++k) {
      buffer[k+offset] = state[k*num_reactor_states + j];
      if(j==num_species+1){
        buffer[k+offset] *= params.ref_temperature_;
      }
      if(j==num_species+2){
        buffer[k+offset] *= params.ref_momentum_;
      }
      if(j==num_species+3){
        buffer[k+offset] *= params.ref_momentum_;
      }
    }

    if(my_pe == npes-1) {
      // Write right (oxidizer) BC data
      if(j==num_species){
        buffer[num_local_points+offset] = params.oxidizer_relative_volume_;
      } else if(j==num_species+1){
        buffer[num_local_points+offset] = params.oxidizer_temperature_*params.ref_temperature_;
      } else if (j==num_species+2) {
        if(params.flame_type_ == 0 || params.flame_type_ == 2) {
          buffer[num_local_points+offset] = 0.0;
        } else if (params.flame_type_ == 1) {
          buffer[num_local_points+offset] = params.G_right_*params.ref_momentum_;
        }
      } else if (j==num_species+3) {
        buffer[num_local_points+offset] = params.P_right_*params.ref_momentum_;
      } else {
        buffer[num_local_points+offset] = params.oxidizer_mass_fractions_[j];
      }
    }

    disp = header_size;
    disp += j*global_size;
    disp += field_index[my_pe]*sizeof(double);
    MPI_File_write_at(output_file, disp, &buffer[0], field_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
  }
  int vars_written = num_reactor_states;

  if(dump_mole_fractions) {
    std::vector<double> inlet_mole_fractions(num_species);
    std::vector<double> fuel_mole_fractions(num_species);
    std::vector<double> oxidizer_mole_fractions(num_species);
    std::vector<double> state_mole_fractions(num_species*num_local_points);
    params.reactor_->GetMechanism()->getXfromY(&params.fuel_mass_fractions_[0],&fuel_mole_fractions[0]);
    params.reactor_->GetMechanism()->getXfromY(&params.inlet_mass_fractions_[0],&inlet_mole_fractions[0]);
    params.reactor_->GetMechanism()->getXfromY(&params.oxidizer_mass_fractions_[0],&oxidizer_mole_fractions[0]);
    for (int k=0; k<num_local_points; ++k) {
      params.reactor_->GetMechanism()->getXfromY(&state[k*num_reactor_states],&state_mole_fractions[k*num_species]);
    }

    for(int j=0; j<num_species; ++j) {
      // Write left (fuel) BC data
      int offset = 0;
      if(my_pe == 0) {
        offset+=1;
        if(params.flame_type_ == 0) {
          buffer[0] = fuel_mole_fractions[j];
        } else {
          buffer[0] = inlet_mole_fractions[j];
        }
      }

      // Write interior data
      for (int k=0; k<num_local_points; ++k) {
        buffer[k+offset] = state_mole_fractions[k*num_species + j];
      }

      if(my_pe == npes-1) {
        // Write right (oxidizer) BC data
        buffer[num_local_points+offset] = oxidizer_mole_fractions[j];
      }

      disp = header_size;
      disp += (vars_written+j)*global_size;
      disp += field_index[my_pe]*sizeof(double);
      MPI_File_write_at(output_file, disp, &buffer[0], field_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }
    vars_written += num_species;
  }

  if(dump_step_rates) {
    std::vector<double> creation_rates(num_steps);
    std::vector<double> destruction_rates(num_steps);
    std::vector<double> concentrations(num_species);

    std::vector<double> inlet_step_rates(num_steps);
    std::vector<double> inlet_net_rates(num_species);
    double relative_volume = params.inlet_relative_volume_;
    double density         = 1.0/relative_volume;
    double temperature     = params.fuel_temperature_*params.ref_temperature_;
    params.reactor_->GetMechanism()->getCfromVY(relative_volume,&params.inlet_mass_fractions_[0],&concentrations[0]);
    params.reactor_->GetMechanism()->getReactionRates(temperature,&concentrations[0],
		                                      &inlet_net_rates[0],
		                                      &creation_rates[0],
		                                      &destruction_rates[0],
		                                      &inlet_step_rates[0]);

    std::vector<double> fuel_step_rates(num_steps);
    std::vector<double> fuel_net_rates(num_species);
    relative_volume = params.fuel_relative_volume_;
    density         = 1.0/relative_volume;
    temperature     = params.fuel_temperature_*params.ref_temperature_;
    params.reactor_->GetMechanism()->getCfromVY(relative_volume,&params.fuel_mass_fractions_[0],&concentrations[0]);
    params.reactor_->GetMechanism()->getReactionRates(temperature,&concentrations[0],
		                                      &fuel_net_rates[0],
		                                      &creation_rates[0],
		                                      &destruction_rates[0],
		                                      &fuel_step_rates[0]);

    std::vector<double> oxidizer_step_rates(num_steps);
    std::vector<double> oxidizer_net_rates(num_species);
    relative_volume = params.oxidizer_relative_volume_;
    density         = 1.0/relative_volume;
    temperature     = params.oxidizer_temperature_*params.ref_temperature_;
    params.reactor_->GetMechanism()->getCfromVY(relative_volume,&params.oxidizer_mass_fractions_[0],&concentrations[0]);
    params.reactor_->GetMechanism()->getReactionRates(temperature,&concentrations[0],
		                                      &oxidizer_net_rates[0],
		                                      &creation_rates[0],
		                                      &destruction_rates[0],
		                                      &oxidizer_step_rates[0]);

    std::vector<double> state_step_rates(num_steps*num_local_points, 0.0);
    std::vector<double> state_net_rates(num_species*num_local_points, 0.0);
    for (int k=0; k<num_local_points; ++k) {
      relative_volume = state[k*num_reactor_states+num_species];
      density         = 1.0/relative_volume;
      temperature     = state[k*num_reactor_states+num_species+1]*params.ref_temperature_;
      params.reactor_->GetMechanism()->getCfromVY(relative_volume,&state[k*num_reactor_states],&concentrations[0]);
      params.reactor_->GetMechanism()->getReactionRates(temperature,&concentrations[0],
                                                        &state_net_rates[k*num_species],
                                                        &creation_rates[0],
                                                        &destruction_rates[0],
                                                        &state_step_rates[k*num_steps]);
    }

    for(int j=0; j<num_steps; ++j) {
      int offset = 0;
      if(my_pe == 0) {
        offset += 1;
        // Write left (fuel) BC data
        if(params.flame_type_ == 0) {
          buffer[0] = fuel_step_rates[j];
        } else {
          buffer[0] = inlet_step_rates[j];
        }
      }

      // Write interior data
      for (int k=0; k<num_local_points; ++k) {
        buffer[k+offset] = state_step_rates[k*num_steps + j];
      }

      // Write right (oxidizer) BC data
      if(my_pe == npes-1) {
        buffer[num_local_points+offset] = oxidizer_step_rates[j];
      }

      disp = header_size;
      disp += (vars_written+j)*global_size;
      disp += field_index[my_pe]*sizeof(double);
      MPI_File_write_at(output_file, disp, &buffer[0], field_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }
    vars_written += num_steps;

    for(int j=0; j<num_species; ++j) {
      int offset = 0;
      if(my_pe == 0) {
        offset += 1;
        // Write left (fuel) BC data
        if(params.flame_type_ == 0) {
          buffer[0] = fuel_net_rates[j];
        } else {
          buffer[0] = inlet_net_rates[j];
        }
      }

      // Write interior data
      for (int k=0; k<num_local_points; ++k) {
        buffer[k+offset] = state_net_rates[k*num_species + j];
      }


      if(my_pe == npes-1) {
        // Write right (oxidizer) BC data
        buffer[num_local_points+offset] = oxidizer_net_rates[j];
      }
      disp = header_size;
      disp += (vars_written+j)*global_size;
      disp += field_index[my_pe]*sizeof(double);
      MPI_File_write_at(output_file, disp, &buffer[0], field_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }
    vars_written += num_species;
  }

  int offset = 0;
  if(my_pe == 0) {
    offset += 1;
    // Write mass flux
    // Left BC
    buffer[0] = params.mass_flux_fuel_;
  }

  // Interior
  for (int k=0; k<num_local_points; ++k) {
    buffer[k+offset] = params.mass_flux_ext_[k+2];
  }

  if(my_pe == npes-1) {
    // Right BC
    buffer[num_local_points+offset] = params.mass_flux_oxidizer_;
  }
  disp = header_size;
  disp += vars_written*global_size;
  disp += field_index[my_pe]*sizeof(double);
  MPI_File_write_at(output_file, disp, &buffer[0], field_count, MPI_DOUBLE, MPI_STATUS_IGNORE);

  MPI_File_close(&output_file);
}


static double GetMixtureMolecularMass(const int grid_id,
                                      const double t,
                                      const double state[],
                                      const FlameParams &params)
{
  const int num_reactor_states = params.reactor_->GetNumStates();
  const int num_species = params.reactor_->GetNumSpecies();
  double molecular_mass=0.0;

  if(0 <= grid_id && grid_id < (int)params.z_.size()) {

    molecular_mass=0.0;
    for(int j=0; j<num_species; ++j) {
      molecular_mass += state[grid_id*num_reactor_states+j]*
        params.inv_molecular_mass_[j];
    }
    molecular_mass = 1.0/molecular_mass;
  }
  return molecular_mass;
}

static double minSumMassFractions(const double state[],
                                  const FlameParams &params)
{
  const int num_reactor_states = params.reactor_->GetNumStates();
  const int num_species = params.reactor_->GetNumSpecies();
  const int num_local_points = (int)params.num_local_points_;
  double sum_mass_fraction;
  double local_min;
  double global_min;

  local_min = 100;
  for(int k=0; k<num_local_points; ++k) {
    sum_mass_fraction=0.0;
    for(int j=0; j<num_species; ++j) {
      sum_mass_fraction += state[k*num_reactor_states+j];
    }
    if (sum_mass_fraction < local_min) {
      local_min = sum_mass_fraction;
    }
  }

  MPI_Allreduce(&local_min,&global_min,1,PVEC_REAL_MPI_TYPE,MPI_MIN,MPI_COMM_WORLD);
  return global_min;
}
static double maxSumMassFractions(const double state[],
                                  const FlameParams &params)
{
  const int num_reactor_states = params.reactor_->GetNumStates();
  const int num_species = params.reactor_->GetNumSpecies();
  const int num_local_points = (int)params.num_local_points_;
  double sum_mass_fraction;
  double local_max, global_max;

  local_max = -100.0;

  for(int k=0; k<num_local_points; ++k) {
    sum_mass_fraction=0.0;
    for(int j=0; j<num_species; ++j) {
      sum_mass_fraction += state[k*num_reactor_states+j];
    }
    if (sum_mass_fraction > local_max) {
      local_max = sum_mass_fraction;
    }
  }

  MPI_Allreduce(&local_max,&global_max,1,PVEC_REAL_MPI_TYPE,MPI_MAX,MPI_COMM_WORLD);
  return global_max;
}
static double minVelocity(const double state[],
			  const FlameParams &params)
{
  const int num_reactor_states = params.reactor_->GetNumStates();
  const int num_local_points = (int)params.num_local_points_;
  double velocity;
  double local_min;
  double global_min;

  local_min = 100000;
  for(int k=0; k<num_local_points; ++k) {
    velocity = params.mass_flux_[k]*state[(k+1)*num_reactor_states-2];
    if (velocity < local_min) {
      local_min = velocity;
    }
  }

  MPI_Allreduce(&local_min,&global_min,1,PVEC_REAL_MPI_TYPE,MPI_MIN,MPI_COMM_WORLD);
  return global_min;
}
static double maxVelocity(const double state[],
			  const FlameParams &params)
{
  const int num_reactor_states = params.reactor_->GetNumStates();
  const int num_local_points = (int)params.num_local_points_;
  double velocity;
  double local_max;
  double global_max;

  local_max = -100000;
  for(int k=0; k<num_local_points; ++k) {
    velocity = params.mass_flux_[k]*state[(k+1)*num_reactor_states-2];
    if (velocity > local_max) {
      local_max = velocity;
    }
  }

  MPI_Allreduce(&local_max,&global_max,1,PVEC_REAL_MPI_TYPE,MPI_MAX,MPI_COMM_WORLD);
  return global_max;
}
static double GetDensityFromYTP(const int grid_id,
                                const double t,
                                const double state[],
                                const FlameParams &params)
{
  const int num_reactor_states = params.reactor_->GetNumStates();
  double density = 0.0;
  double molecular_mass = GetMixtureMolecularMass(grid_id,
                                                  t,
                                                  state,
                                                  params);
  double temperature;

  if(0 <= grid_id && grid_id < (int)params.z_.size()) {

    temperature = state[(grid_id+1)*num_reactor_states-1]*
      params.ref_temperature_;
    density = params.parser_->pressure()*molecular_mass/
      (params.reactor_->GetGasConstant()*temperature);

  }
  return density;
}

static int GetTrackMaxStateId(const FlameParams &params,
                              std::vector<int> *state_id)
{
  std::string prefix="MassFraction_";
  std::string state_name;
  const int num_reactor_states = params.reactor_->GetNumStates();
  const int num_species = params.reactor_->GetNumSpecies();
  int species_id;
  state_id->clear();
  state_id->push_back(num_species+1); // set the temperature id

  for(size_t j=0; j<params.parser_->track_max().size(); ++j) {
    state_name.clear();
    state_name = prefix + params.parser_->track_max()[j];
    species_id = params.reactor_->GetIdOfState(state_name.c_str());
    if(0 <= species_id && species_id < num_species) {
      state_id->push_back(species_id);
    } else {
      params.logger_->PrintF(
        "# WARNING: In GetTrackMaxStateId,\n"
        "#          input file species given by track_max[%lu] = '%s'\n"
        "#          is not found in the mechanism. Skipping species.\n",
        j,
        params.parser_->track_max()[j].c_str());
      params.logger_->FFlush();
    }
  }
  return (int)state_id->size();
}

static int GetStateMaxima(const std::vector<int> &state_id,
                          const double state[],
                          const FlameParams &params,
                          std::vector<double> *max_value,
                          std::vector<double> *max_position)
{
  const int num_local_points = (int)params.num_local_points_;
  const int num_reactor_states = params.reactor_->GetNumStates();
  const int num_species = params.reactor_->GetNumSpecies();
  double state_max, state_max_position;
  max_value->clear();
  max_position->clear();

  for(size_t j=0; j<state_id.size(); ++j) {

    state_max = FindMaximumParallel(num_local_points,
				    &params.z_[0],
				    1,
				    &state[state_id[j]],
				    num_reactor_states,
				    true, // use quadratic
				    &state_max_position);
    if(state_id[j] == num_species+1) {
      // temperature
      state_max*=params.ref_temperature_;
    }
    max_value->push_back(state_max);
    max_position->push_back(state_max_position);
  }

  return 0;
}

static int GetFuelSpeciesId(const FlameParams &params,
                            std::vector<int> *species_id)
{
  const int num_species = params.reactor_->GetNumSpecies();
  std::map<std::string, double>::const_iterator iter;
  species_id->clear();
  for(iter=params.parser_->inlet_fuel_comp().begin();
      iter != params.parser_->inlet_fuel_comp().end();
      ++iter) {

    std::string species_name=iter->first;
    std::string state_name = std::string("MassFraction_")+species_name;
    double mole_fraction = iter->second;
    int id = params.reactor_->GetIdOfState(state_name.c_str());
    if(0 <= id && id < num_species && mole_fraction > 0.0) {
      species_id->push_back(id);
    }
  }
  return (int)species_id->size();
}

static double InletFuelFraction(const FlameParams &params)
{
  double fuel_sum = 0.0;
  std::vector<int> fuel_species_id;
  GetFuelSpeciesId(params,&fuel_species_id);
  for(size_t j=0; j<fuel_species_id.size(); ++j) {
    fuel_sum += params.fuel_mass_fractions_[fuel_species_id[j]];
  }
  return fuel_sum;
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
*/

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}

static double StoichiometricMixtureFraction(const FlameParams &params)
{
  double fuel_sum = 0.0;
  std::vector<int> fuel_species_id;
  GetFuelSpeciesId(params,&fuel_species_id);
  for(size_t j=0; j<fuel_species_id.size(); ++j) {
    fuel_sum += params.stoichiometric_mass_fractions_[fuel_species_id[j]];
  }
  return fuel_sum;
}

static int SootOutput(const FlameParams &params,
                       const double state[])
{
  const int num_species = params.reactor_->GetNumSpecies();
  const int num_states = params.reactor_->GetNumStates();
  const int num_local_points = params.num_local_points_;
  std::vector<string> state_name;
  std::vector<double> diameter;
  std::vector<double> mole_fraction;
  std::vector<double> number_density;
  std::vector<double> soot_vol_frac;
  int num_tracked_species;
  string prefix = "MassFraction_";
  string line;
  ifstream PAH_file;
  double mixture_molecular_mass;
  int ind;

  // Get list of tracked species
  PAH_file.open("PAH_file");
  if(params.my_pe_ == 0) {
    if(!PAH_file.is_open()) {
      printf("# Error opening PAH_file\n");
      return 1;
    }
  }

  while(getline(PAH_file,line)) {
    istringstream ss(line);
    std::string word;
    int count = 0;
    while(ss >> word) {
      if(count == 0) {
        word = prefix + word;
        state_name.push_back(word);
      }
      if(count == 1) {
        diameter.push_back(std::stod(word));
      }
      count++;
    }
  }
  num_tracked_species = state_name.size();

  mole_fraction.assign(num_tracked_species*num_local_points, 0.0);
  number_density.assign(num_tracked_species*num_local_points, 0.0);
  soot_vol_frac.assign(num_local_points, 0.0);

  // Compute mole fractions and number densities of tracked species
  for(int j=0; j<num_local_points; j++) {

    mixture_molecular_mass = 0.0;
    for(int k=0; k<num_species; k++) {
      mixture_molecular_mass += state[j*num_states + k]*params.inv_molecular_mass_[k];
    }
    mixture_molecular_mass = 1.0/mixture_molecular_mass;

    soot_vol_frac[j] = 0.0;
    for(int i=0; i<num_tracked_species; i++) {
      if((ind = params.reactor_->GetIdOfState(state_name[i].c_str()) )  != -1) {

        mole_fraction[j*num_tracked_species+i] = state[j*num_states + ind]*mixture_molecular_mass*
          params.inv_molecular_mass_[ind];
        // n = N_A*conc = N_A*rho*Yi/Wi
        number_density[j*num_tracked_species+i] = params.mechanism_->getAvogadroNumber()/
          state[j*num_states+num_species]*state[j*num_states + ind]*
          params.inv_molecular_mass_[ind];

        soot_vol_frac[j] += number_density[j*num_tracked_species+i]*
          pow(diameter[i],3.0)*3.1416/6;
      } // if tracked_species
    } // loop tracked_species

    // Compute soot volume fraction

  } // loop grid

  // Output
  int npes = params.npes_;
  int my_pe = params.my_pe_;
  MPI_Comm comm = params.comm_;
  std::vector<double> mole_fraction_all, number_density_all, soot_vol_frac_all;
  mole_fraction_all.assign(num_local_points*npes*num_tracked_species, 0.0);
  number_density_all.assign(num_local_points*npes*num_tracked_species, 0.0);
  soot_vol_frac_all.assign(num_local_points*npes, 0.0);

  long int dsize = num_local_points*num_tracked_species;
  int nodeDest = 0;
  MPI_Gather(&mole_fraction[0], dsize, PVEC_REAL_MPI_TYPE,
             &mole_fraction_all[0], dsize, PVEC_REAL_MPI_TYPE, nodeDest, comm);

  MPI_Gather(&number_density[0], dsize, PVEC_REAL_MPI_TYPE,
             &number_density_all[0], dsize, PVEC_REAL_MPI_TYPE, nodeDest, comm);

  dsize = num_local_points;
  MPI_Gather(&soot_vol_frac[0], dsize, PVEC_REAL_MPI_TYPE,
             &soot_vol_frac_all[0], dsize, PVEC_REAL_MPI_TYPE, nodeDest, comm);


  if(my_pe == 0) {
    FILE * dimerFile;
    dimerFile = fopen("soot_output","w");
    fprintf(dimerFile,"# x (m)   soot_VF  ");
    for(int i=0; i<num_tracked_species; i++) {
      line = "X_" + state_name[i].substr(state_name[i].find("_") + 1);
      fprintf(dimerFile,"%s  ",line.c_str());
      line = "n_" + state_name[i].substr(state_name[i].find("_") + 1);
      fprintf(dimerFile,"%s  ",line.c_str());
    }
    fprintf(dimerFile,"\n");
    for(int j=0; j<num_local_points*npes; j++) {
      fprintf(dimerFile,"%14.7e  %14.7e  ",params.z_[j], soot_vol_frac_all[j]);
      for(int i=0; i<num_tracked_species; i++) {
        fprintf(dimerFile,"%14.7e   %14.7e  ",
                mole_fraction_all[j*num_tracked_species+i],
                number_density_all[j*num_tracked_species+i]);
      }
      fprintf(dimerFile,"\n");
    }
    fclose(dimerFile);
  }

  return 0;
}
