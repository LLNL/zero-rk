#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <string>
#include <map>

#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
#ifdef ZERORK_MPI
#include <nvector/nvector_parallel.h> // parallel N_Vector types, fcts., and macros
#else
#include <nvector/nvector_serial.h>
#endif

#ifdef SUNDIALS2
#include <cvode/cvode_spgmr.h>
#elif SUNDIALS3
#include <cvode/cvode_spils.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#elif SUNDIALS4
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif

#include <utilities/file_utilities.h>

#ifdef ZERORK_MPI
#include <mpi.h>
#endif

#include "flame_params.h"
#include "cvode_functions.h"
#include "set_initial_conditions.h"

using zerork::getHighResolutionTime;

const bool RUN_DEBUG=false;
const int NUM_STDOUT_PARAMS = 19; // number of non species parameters to write
                                  // to standard out

static int check_flag(void *flagvalue, const char *funcname, int opt);

static double FindMaximum(const size_t num_points,
                          const double x[],
                          const size_t x_stride,
                          const double f[],
                          const size_t f_stride,
                          const bool use_quadratic,
                          double *x_at_max);
#ifdef ZERORK_MPI
static void WriteFieldParallel(double t,
			       const double state[],
			       const FlameParams &params);
#endif
static void WriteFieldSerial(double t,
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
#ifdef ZERORK_MPI
  MPI_Init(&argc, &argv);
#endif

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
  FlameParams flame_params(argv[1]);

  // Initialize state variables vector and pointer
  N_Vector flame_state;
  flame_state = NULL;
  double *flame_state_ptr;
  flame_state_ptr = NULL;

  // Get constants from flame params
  const int num_grid_points = flame_params.z_.size();
  const int num_local_points = flame_params.num_local_points_;
  const int num_states = flame_params.reactor_->GetNumStates();
  long int num_local_states = num_local_points*num_states;
  long int total_states = num_grid_points*num_states;
  double max_time = flame_params.parser_->max_time();
  const double step_dt = flame_params.parser_->stats_dt();
  const double ref_temperature = flame_params.ref_temperature_;
  const double dz = flame_params.dz_[(int)num_grid_points/2]; //should be min(dz_)
  int my_pe = flame_params.my_pe_;

  // Declare variables
  int flag=0;
  int num_prints = 0;
  double current_time = 0.0;
  double time_offset;
  double next_time;
  double max_temperature_jump, z_max_temperature_jump;
  double inlet_density, inlet_molecular_mass;
  double inlet_fuel_fraction;
  double min_sum_mass_fraction, max_sum_mass_fraction;
  double min_velocity, max_velocity;

  // CVode integrator Stats
  long int nsteps, nfevals, nlinsetups, netfails;
  int qlast, qcur;
  double hinused, hlast, hcur, tcur;
  zerork::utilities::Logger field_file(flame_params.parser_->field_file(),
                                        false,
                                        false);
  std::vector<double> temperature_jump;
  std::vector<int> track_max_state_id;
  std::vector<double> state_maxima, state_maxima_positions;

  // allocate CVode data structures
#ifdef ZERORK_MPI
  flame_state     = N_VNew_Parallel(MPI_COMM_WORLD, num_local_states, total_states);
  flame_state_ptr = NV_DATA_P(flame_state);
#else
  flame_state     = N_VNew_Serial(total_states);
  flame_state_ptr = NV_DATA_S(flame_state);
#endif

  // Initialize state vector
  time_offset = 0.0;
  SetInitialCompositionAndWallTemp(flame_params, flame_state_ptr, &time_offset);
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

  flag = CVodeSetLinearSolver(cvode_ptr, LS, NULL);
  if(check_flag(&flag, "CVodeSetLinearSolver", 1)) exit(-1);

  flag = SUNSPGMRSetGSType(LS, MODIFIED_GS);
  if(check_flag(&flag, "SUNSPGMRSetGSType", 1)) exit(-1);

  /* Set the preconditioner setup and solve functions */
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
#ifdef ZERORK_MPI
    WriteFieldParallel(time_offset,
		       &flame_state_ptr[0],
		       flame_params);
#else
    WriteFieldSerial(time_offset,
                     &flame_state_ptr[0],
                     flame_params);
#endif
  }

  // Report simulation info
  if(my_pe == 0) {
    inlet_molecular_mass = GetMixtureMolecularMass(
	   0,
	   0.0,
	   &flame_params.inlet_mass_fractions_[0],
	   flame_params);
    inlet_density = flame_params.parser_->pressure()*inlet_molecular_mass/
      (flame_params.inlet_temperature_*flame_params.ref_temperature_*
       flame_params.reactor_->GetGasConstant());
    inlet_fuel_fraction = InletFuelFraction(flame_params);
    printf("# Number of states     : %d\n",num_states);
    printf("# Number of grid points: %d\n", num_grid_points);
    printf("# Inlet BC fuel fraction [kg fuel/kg inflow]: %.18g\n",
	   inlet_fuel_fraction);
    printf("# Inlet BC temperature  [K]: %.18g\n",
	   flame_params.inlet_temperature_*flame_params.ref_temperature_);
    printf("# Inlet BC density [kg/m^3]: %.18g\n", inlet_density);
    printf("# Inlet BC velocity   [m/s]: %.18g\n", flame_params.mass_flux_[0]/inlet_density);
    printf("# Initial upstream temperature  (next to inlet)     [K]: %.18g\n",
	   flame_state_ptr[num_states-1]*ref_temperature);
    printf("# Initial upstream density      (next to inlet)[kg/m^3]: %.18g\n",
	   1.0/flame_state_ptr[num_states-2]);
    printf("# Initial upstream avg velocity (next to inlet)   [m/s]: %.18g\n",
	   flame_params.mass_flux_[0]*flame_state_ptr[num_states-2]);
    printf("#------------------------------------------------------------------------------\n");
    printf("# Column  1: [#] time steps printed\n");
    printf("# Column  2: [s] ODE system time\n");
    printf("# Column  3: [m/s] Flame speed\n");
    printf("# Column  4: [m] Flame thickness defined as T_burnt-T_unburnt\n"
	   "#            divided by maximum temperature gradient\n");
    printf("# Column  5: [m] Flame thickness defined as (K/rhoCp)/SL\n");
    printf("# Column  6: [m] location of the maximum temperature \n");
    printf("# Column  7: [#] number of steps taken by CVode\n");
    printf("# Column  8: [#] number of calls to the user's time derivative (RHS) function\n");
    printf("# Column  9: [#] number of calls to the linear solver setup function\n");
    printf("# Column 10: [#] number of error test failures\n");
    printf("# Column 11: [#] current method order of time integrator\n");
    printf("# Column 12: [s] current size of internal time step\n");
    printf("# Column 13: [s] characteristic time step for Courant number control dx/max(|velocity|) \n");
    printf("# Column 14: [s] characteristic time step for diffusion stability control 0.5*dx*dx/max(thermal diffusivity)\n");
    printf("# Column 15: Mass change normalized by inlet mass flux\n");
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

      // Compute T-Twall
      for(int j=0; j<num_local_points; ++j) {
	int jglobal = j + my_pe*num_local_points;
        temperature_jump[j] =
          (flame_state_ptr[(j+1)*num_states-1]-
           flame_params.wall_temperature_[jglobal])*flame_params.ref_temperature_;

      }
      // Find maximum value over domain
      max_temperature_jump = FindMaximum(num_local_points,
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

      // Print to screen/log
      if(my_pe == 0) { //only root prints
	printf("%7d  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %6d  %6d  %6d  %6d  %6d  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
	       num_prints,
	       current_time+time_offset,
	       flame_params.flame_speed_,
	       flame_params.flame_thickness_,
	       flame_params.flame_thickness_alpha_,
	       z_max_temperature_jump,
	       (int)nsteps,
	       (int)nfevals,
	       (int)nlinsetups,
	       (int)netfails,
	       qcur,
	       hcur,
	       dz/flame_params.max_velocity_,
	       0.5*dz*dz/flame_params.max_thermal_diffusivity_,
	       flame_params.mass_change_,
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
#ifdef ZERORK_MPI
	WriteFieldParallel(current_time+time_offset,
			   &flame_state_ptr[0],
			   flame_params);
#else
        WriteFieldSerial(current_time+time_offset,
                         &flame_state_ptr[0],
                         flame_params);
#endif
      }


    } else {
      printf("# ERROR: In time marching loop,\n");
      printf("#        CVode returned flag = %d.\n", flag);
      printf("#        Halting integration.\n");
      break;
    }

  }
  loop_time = getHighResolutionTime() - clock_time;
#ifdef ZERORK_MPI
  N_VDestroy_Parallel(flame_state);
#else
  N_VDestroy_Serial(flame_state);
#endif
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

#ifdef ZERORK_MPI
  MPI_Finalize();
#endif
  return 0;
}

// Find maximum and its location in parallel
static double FindMaximum(const size_t num_points,
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

#ifdef ZERORK_MPI
  // Compute global maximum
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  in.index += myrank*num_points;

  MPI_Allreduce(&in,&out,1,MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
#else
  out.index = in.index;
  out.value = in.value;
#endif

  *x_at_max = x[out.index];
  return out.value;
}

// Write field to binary file for restart/postprocessing
#ifdef ZERORK_MPI
static void WriteFieldParallel(double t,
			       const double state[],
			       const FlameParams &params)
{
  int num_grid_points = (int)params.z_.size();
  int num_grid_points_ext = (int)params.z_.size() + 1; // left BC
  int num_local_points = (int)params.num_local_points_;
  int num_reactor_states = params.reactor_->GetNumStates();
  int my_pe, npes;
  char filename[32], *basename;

  int disp;
  std::vector<double> buffer;
  buffer.assign(num_local_points, 0.0);

  MPI_File output_file;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_pe);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  MPI_Datatype localarray;
  int order = MPI_ORDER_C;
  int gsize[1] = {num_grid_points};
  int lsize[1] = {num_local_points};
  int start[1] = {my_pe*num_local_points};

  MPI_Type_create_subarray(1,gsize,lsize,start,
			   order,MPI_DOUBLE,&localarray);
  MPI_Type_commit(&localarray);

  basename = "data_";
  sprintf(filename, "%s%f", basename, t);

  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_RDWR,
		MPI_INFO_NULL, &output_file);

  // Write header
  if(my_pe == 0) {
    MPI_File_write(output_file, &num_grid_points_ext, 1, MPI_INT, MPI_STATUS_IGNORE);//num points
    MPI_File_write(output_file, &num_reactor_states, 1, MPI_INT, MPI_STATUS_IGNORE);//num variables
    MPI_File_write(output_file, &t, 1, MPI_DOUBLE, MPI_STATUS_IGNORE); //time
    for(int j=0; j<num_reactor_states; ++j) {
      std::string state_name = params.reactor_->GetNameOfStateId(j);
      char buf[64];
      strcpy(buf, state_name.c_str());
      MPI_File_write(output_file, &buf, 64, MPI_CHAR, MPI_STATUS_IGNORE); //state name
    }
  }

  // Write data for each variable
  for(int j=0; j<num_reactor_states; ++j) {
    // Left BC
    if(j==num_reactor_states-1) {
      buffer[0] = params.inlet_temperature_*params.ref_temperature_;
    } else if(j==num_reactor_states-2) {
      buffer[0] = params.inlet_relative_volume_;
    } else {
      buffer[0] = params.inlet_mass_fractions_[j];
    }
    disp = 2*sizeof(int) + sizeof(double) + num_reactor_states*sizeof(char)*64
      + j*(npes*num_local_points+1)*sizeof(double);
    MPI_File_set_view(output_file, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    MPI_File_write_all(output_file, &buffer[0], 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

    // Write interior data
    for (int k=0; k<num_local_points; ++k) {
      buffer[k] = state[k*num_reactor_states + j];
      if(j==num_reactor_states-1){
	buffer[k] *= params.ref_temperature_;
      }
    }
    disp = 2*sizeof(int) + sizeof(double) + num_reactor_states*sizeof(char)*64
      + sizeof(double)
      + j*(npes*num_local_points+1)*sizeof(double);
    MPI_File_set_view(output_file, disp, MPI_DOUBLE, localarray, "native", MPI_INFO_NULL);
    MPI_File_write_all(output_file, &buffer[0], num_local_points, MPI_DOUBLE, MPI_STATUS_IGNORE);
  }
  // Write mass flux
  // Left BC
  buffer[0] = params.mass_flux_[0];
  disp = 2*sizeof(int) + sizeof(double) + num_reactor_states*sizeof(char)*64
    + num_reactor_states*(npes*num_local_points+1)*sizeof(double);
  MPI_File_set_view(output_file, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
  MPI_File_write_all(output_file, &buffer[0], 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

  // Interior
  for (int k=0; k<num_local_points; ++k) {
    buffer[k] = params.mass_flux_[k];
  }
  disp = 2*sizeof(int) + sizeof(double) + num_reactor_states*sizeof(char)*64
    + sizeof(double)
    + num_reactor_states*(npes*num_local_points+1)*sizeof(double);
  MPI_File_set_view(output_file, disp, MPI_DOUBLE, localarray, "native", MPI_INFO_NULL);
  MPI_File_write_all(output_file, &buffer[0], num_local_points, MPI_DOUBLE, MPI_STATUS_IGNORE);

  MPI_File_close(&output_file);

  MPI_Type_free(&localarray);
}
#endif

static void WriteFieldSerial(double t,
                             const double state[],
                             const FlameParams &params)
{
  int num_grid_points = (int)params.z_.size();
  int num_grid_points_ext = (int)params.z_.size() + 1; // left BC
  int num_local_points = (int)params.num_local_points_;
  int num_reactor_states = params.reactor_->GetNumStates();
  char filename[32], *basename;

  int disp;
  std::vector<double> buffer;
  buffer.assign(num_local_points, 0.0);

  ofstream output_file;

  basename = "data_";
  sprintf(filename, "%s%f", basename, t);

  output_file.open(filename, ios::binary | ios::out);

  // Write header
  output_file.write((char*)&num_grid_points_ext, sizeof(int));
  output_file.write((char*)&num_reactor_states, sizeof(int));
  output_file.write((char*)&t, sizeof(double));
  for(int j=0; j<num_reactor_states; j++) {
    std::string state_name = params.reactor_->GetNameOfStateId(j);
    char buf[64];
    strcpy(buf, state_name.c_str());
    output_file.write((char*)&buf, 64*sizeof(char));
  }

  // Write data for each variable
  for(int j=0; j<num_reactor_states; ++j) {
    // Left BC
    if(j==num_reactor_states-1) {
      buffer[0] = params.inlet_temperature_*params.ref_temperature_;
    } else if(j==num_reactor_states-2) {
      buffer[0] = params.inlet_relative_volume_;
    } else {
      buffer[0] = params.inlet_mass_fractions_[j];
    }
    output_file.write((char*)&buffer[0], sizeof(double));

    // Interior data
    for (int k=0; k<num_local_points; ++k) {
      buffer[k] = state[k*num_reactor_states + j];
      if(j==num_reactor_states-1){
	buffer[k] *= params.ref_temperature_;
      }
    }
    output_file.write((char*)&buffer[0], num_local_points*sizeof(double));
  }

  // Write mass flux
  // Left BC
  buffer[0] = params.mass_flux_[0];
  output_file.write((char*)&buffer[0], sizeof(double));
  // Interior
  for (int k=0; k<num_local_points; ++k) {
    buffer[k] = params.mass_flux_[k];
  }
  output_file.write((char*)&buffer[0], num_local_points*sizeof(double));
  output_file.close();

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
#ifdef ZERORK_MPI
  MPI_Allreduce(&local_min,&global_min,1,PVEC_REAL_MPI_TYPE,MPI_MIN,MPI_COMM_WORLD);
#else
  global_min = local_min;
#endif
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
#ifdef ZERORK_MPI
  MPI_Allreduce(&local_max,&global_max,1,PVEC_REAL_MPI_TYPE,MPI_MAX,MPI_COMM_WORLD);
#else
  global_max = local_max;
#endif
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
#ifdef ZERORK_MPI
  MPI_Allreduce(&local_min,&global_min,1,PVEC_REAL_MPI_TYPE,MPI_MIN,MPI_COMM_WORLD);
#else
  global_min = local_min;
#endif
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
#ifdef ZERORK_MPI
  MPI_Allreduce(&local_max,&global_max,1,PVEC_REAL_MPI_TYPE,MPI_MAX,MPI_COMM_WORLD);
#else
  global_max = local_max;
#endif
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
  state_id->push_back(num_reactor_states-1); // set the temperature id

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

  double state_max, state_max_position;
  max_value->clear();
  max_position->clear();

  for(size_t j=0; j<state_id.size(); ++j) {

    state_max = FindMaximum(num_local_points,
                            &params.z_[0],
                            1,
                            &state[state_id[j]],
                            num_reactor_states,
                            true, // use quadratic
                            &state_max_position);

    if(state_id[j] == num_reactor_states-1) {
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
    fuel_sum += params.inlet_mass_fractions_[fuel_species_id[j]];
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
