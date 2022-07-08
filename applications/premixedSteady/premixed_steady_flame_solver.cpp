#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <string>
#include <map>

#include <utilities/file_utilities.h>

#include <mpi.h>

#include "flame_params.h"
#include "kinsol_functions.h"
#include "set_initial_conditions.h"


#include <nvector/nvector_parallel.h> // serial N_Vector types, fcts., and macros
#include <kinsol/kinsol.h>

#ifdef SUNDIALS2
#include <kinsol/kinsol_spgmr.h>
#elif SUNDIALS3
#include <kinsol/kinsol_spils.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#elif SUNDIALS4
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif
//#include <kinsol/kinsol_bbdpre.h>
//#include <sundials/sundials_dense.h>
//#include <sundials/sundials_types.h>
//#include <sundials/sundials_math.h>

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

int main(int argc, char *argv[])
{
  double clock_time = getHighResolutionTime();
  double setup_time, loop_time, sensanal_time;

  if(argc < 2) {
    printf("# ERROR: Incorrect command line usage.\n");
    printf("#        use %s <input parameters>\n",argv[0]);
    exit(-1);
  }

  // MPI
  MPI_Comm comm;
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  
  // KINSOL memory pointer and linear solver
  void *kinsol_ptr = NULL;
#if defined SUNDIALS3 || defined SUNDIALS4
  SUNLinearSolver LS;
  LS = NULL;
#endif

  // Initialize flame params
  FlameParams flame_params(argv[1],comm);

  // Initialize state variables/scalers vectors and pointers
  N_Vector flame_state, scaler, constraints;
  flame_state = NULL;
  scaler = constraints = NULL;
  double *scaler_ptr;
  scaler_ptr = NULL;
  double *flame_state_ptr;
  flame_state_ptr = NULL;

  //Declare variables
  int Nlocal;

  int maxl, maxlrst;
  int maxiter, mset;
  int flag = 0;
  const int num_grid_points = flame_params.z_.size();
  const int num_local_points = flame_params.num_local_points_;
  const int num_states = flame_params.reactor_->GetNumStates();
  const int num_steps = flame_params.reactor_->GetNumSteps();
  const int num_reactions = flame_params.reactor_->GetNumReactions();
  long int num_local_states = num_local_points*num_states;
  long int total_states = num_grid_points*num_states;
  const double ref_temperature = flame_params.ref_temperature_;
  const double dz = flame_params.dz_[(int)num_grid_points/2]; //should be min(dz_)
  int num_prints = 0;
  double current_time = 0.0;
  double time_offset;
  double max_temperature_jump, z_max_temperature_jump;
  double inlet_density, inlet_molecular_mass;
  double inlet_fuel_fraction;
  double min_sum_mass_fraction, max_sum_mass_fraction;
  double min_velocity, max_velocity;

  Nlocal = num_local_states;

  // KINSol integrator Stats
  long int nsteps, nfevals, nliniters, njacsetups, njacsolves;
  zerork::utilities::Logger field_file(flame_params.parser_->field_file(),
                               false,
                               false);
  std::vector<double> temperature_jump;
  std::vector<int> track_max_state_id;
  std::vector<double> state_maxima, state_maxima_positions;

  int my_pe = flame_params.my_pe_;

  // allocate KINSOL data structures
  flame_state          = N_VNew_Parallel(comm, Nlocal, total_states);
  flame_state_ptr      = NV_DATA_P(flame_state);
  scaler          = N_VNew_Parallel(comm, Nlocal, total_states);
  scaler_ptr      = NV_DATA_P(scaler);
  constraints     = N_VNew_Parallel(comm, Nlocal, total_states);
  N_VConst(0.0, constraints); //0.0 no constraints, 1.0 u_i >= 0.0, 2.0 u_i > 0.0

  // Initialize state vector
  time_offset = 0.0;
  SetInitialCompositionAndWallTemp(flame_params, flame_state_ptr, &time_offset);

  std::vector<int> track_max_all;
  std::vector<double> state_max_all, state_max_pos_all;
  track_max_all.assign(num_states, 0);
  for (int k=0; k<num_states; ++k){
    track_max_all[k] = k;
  }
  GetStateMaxima(track_max_all,
		 flame_state_ptr,
		 flame_params,
		 &state_max_all,
		 &state_max_pos_all);

  for (int j=0; j<num_local_points; ++j) {
    for (int k=0; k<num_states; ++k){
      scaler_ptr[j*num_states + k] = 1.0;
    }
  }

  //----------------------------------------------------------------------------
  // Setup KINSOL solver
  // Create KINSOL pointer
  kinsol_ptr = KINCreate();

  // Initialize KINSOL module with RHS function and state vector
  flag = KINInit(kinsol_ptr, ConstPressureFlame, flame_state);

  // Set function to handle errors and exit cleanly
  flag = KINSetErrHandlerFn(kinsol_ptr, ErrorFunction, &flame_params);

  // Set user data
  flag = KINSetUserData(kinsol_ptr, &flame_params);

  // Set constraints
  flag = KINSetConstraints(kinsol_ptr, constraints);
  N_VDestroy_Parallel(constraints);

  // Set tolerances
  // RHS(y) < fnormtol
  flag = KINSetFuncNormTol(kinsol_ptr, flame_params.parser_->rel_tol());
  // Step tolerance: dy < steptol
  flag = KINSetScaledStepTol(kinsol_ptr, flame_params.parser_->abs_tol());

  // Setup KINSol
  // max Krylov dimension
  maxl = 1000; //TO DO: set it from input file
  maxlrst = 0;

#ifdef SUNDIALS2
    // Initialize Linear Solver
  flag = KINSpgmr(kinsol_ptr, maxl);
  flag = KINSpilsSetMaxRestarts(kinsol_ptr, maxlrst);
  // Set preconditioner
  if (flame_params.integrator_type_ == 2) {
    flag = KINSpilsSetPreconditioner(kinsol_ptr, ReactorBBDSetup, ReactorBBDSolve);
  } else if(flame_params.integrator_type_ == 3){
    flag = KINSpilsSetPreconditioner(kinsol_ptr, ReactorAFSetup, ReactorAFSolve);
  } else {
    printf("integrator_type == %d not currently supported\n",
           flame_params.integrator_type_);
    flag = -1;
  }
#elif SUNDIALS3
  // Initialize Linear Solver
  LS = SUNSPGMR(flame_state, PREC_RIGHT, maxl);
  flag = KINSpilsSetLinearSolver(kinsol_ptr, LS);
  flag = SUNSPGMRSetMaxRestarts(LS, maxlrst);

  // Set preconditioner
  if (flame_params.integrator_type_ == 2) {
    flag = KINSpilsSetPreconditioner(kinsol_ptr, ReactorBBDSetup, ReactorBBDSolve);
  } else if(flame_params.integrator_type_ == 3){
    flag = KINSpilsSetPreconditioner(kinsol_ptr, ReactorAFSetup, ReactorAFSolve);
  } else {
    printf("integrator_type == %d not currently supported\n",
           flame_params.integrator_type_);
    flag = -1;
  }
#elif SUNDIALS4
  // Initialize Linear Solver
  LS = SUNLinSol_SPGMR(flame_state, PREC_RIGHT, maxl);
  flag = KINSetLinearSolver(kinsol_ptr, LS, NULL);
  flag = SUNLinSol_SPGMRSetMaxRestarts(LS, maxlrst);

  // Set preconditioner
  if (flame_params.integrator_type_ == 2) {
    flag = KINSetPreconditioner(kinsol_ptr, ReactorBBDSetup, ReactorBBDSolve);
  } else if(flame_params.integrator_type_ == 3){
    flag = KINSetPreconditioner(kinsol_ptr, ReactorAFSetup, ReactorAFSolve);
  } else {
    printf("integrator_type == %d not currently supported\n",
           flame_params.integrator_type_);
    flag = -1;
  }

#endif

  maxiter = 100;
  flag = KINSetNumMaxIters(kinsol_ptr, maxiter);

  //0 for default, 1 for exact Newton, > 1 for modified Newton
  mset = flame_params.parser_->max_subiter();
  flag = KINSetMaxSetupCalls(kinsol_ptr, mset);

  flag = KINSetPrintLevel(kinsol_ptr, 1);
  // 0 for no info, 1 for scaled l2 norm, 3 for additional linear solver info

  //----------------------------------------------------------------------------

  temperature_jump.assign(num_local_points,0.0);
  GetTrackMaxStateId(flame_params, &track_max_state_id);

  if (time_offset == 0.0) {
    WriteFieldParallel(time_offset,
		       &flame_state_ptr[0],
		       flame_params);
  }

  // Report simulation info to screen/logfile
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
    printf("# Inlet BC velocity   [m/s]: %.18g\n", flame_state_ptr[num_states-2]/inlet_density);
    printf("# Initial upstream temperature  (next to inlet)     [K]: %.18g\n",
	   flame_state_ptr[num_states-1]*ref_temperature);
    printf("# Initial upstream density      (next to inlet)[kg/m^3]: %.18g\n",
	   inlet_density);
    printf("# Initial upstream avg velocity (next to inlet)   [m/s]: %.18g\n",
	   flame_state_ptr[num_states-2]/inlet_density);
    printf("#------------------------------------------------------------------------------\n");
    printf("# Column  1: [#] time steps printed\n");
    printf("# Column  2: [s] System time\n");
    printf("# Column  3: [m/s] Flame speed\n");
    printf("# Column  4: [m] Flame thickness defined as T_burnt-T_unburnt\n"
	   "#            divided by maximum temperature gradient\n");
    printf("# Column  5: [m] Flame thickness defined as (K/rhoCp)/SL\n");
    printf("# Column  6: [m] location of the maximum temperature jump\n");
    printf("# Column  7: [#] number of steps taken by KINSOL\n");
    printf("# Column  8: [#] number of calls to the user's (RHS) function\n");
    printf("# Column  9: [#] number of calls to the linear solver setup function\n");
    printf("# Column 10: [#] number of linear iterations\n");
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

  // Get time
  setup_time = getHighResolutionTime() - clock_time;
  clock_time = getHighResolutionTime();

  double fnorm, stepnorm;

    // Pseudo-unsteady
  if(flame_params.pseudo_unsteady_) {
    flag = KINSetPrintLevel(kinsol_ptr, 0);
    flag = KINSetNumMaxIters(kinsol_ptr, 80);
    if(my_pe==0) {
      printf("# begin pseudo time-stepping\n");
      printf("# time      timestep     SL     nfevals  cputime\n");
    }
    double pseudo_time = 0.0;
    double kinstart_time, kinend_time;
    flame_params.dt_ = flame_params.parser_->pseudo_unsteady_dt()*0.5; //Half the timestep for first iteration

    while(pseudo_time < 0.05) {
      for(int j=0; j<num_local_states; j++) {
        flame_params.y_old_[j] = flame_state_ptr[j];
      }

      // Solve system
      kinstart_time = getHighResolutionTime();
      flag = KINSol(kinsol_ptr,
                    flame_state,
                    KIN_NONE,
                    scaler,
                    scaler);
      kinend_time = getHighResolutionTime();

      KINGetNumFuncEvals(kinsol_ptr,&nfevals);
      KINGetFuncNorm(kinsol_ptr, &fnorm);

      if((flag==0 || flag==1) && fnorm != 0.0){
        // Success
        pseudo_time += flame_params.dt_;
        if(my_pe==0) {printf("%5.3e   %5.3e   %5.3e  %d   %5.3e\n",
                             pseudo_time,
                             flame_params.dt_,
                             flame_params.flame_speed_,
                             nfevals,
                             kinend_time-kinstart_time);}

        // Increase time step if it's converging quickly
        if(nfevals <= 5 || pseudo_time == flame_params.dt_) {
          flame_params.dt_ *= 3.0;
        } else if (nfevals <= 10) {
          flame_params.dt_ *= 1.3;
        }
      } else {
        // Failure, reset state and decrease timestep
        for(int j=0; j<num_local_states; j++) {
          flame_state_ptr[j] = flame_params.y_old_[j];
        }
        flame_params.dt_ *= 0.5;
        if(flame_params.dt_ < 1.0e-5){
          if(my_pe==0) {printf("# Error, pseudo-unsteady solver is not converging\n");}
          exit(-1);
        }
        if(my_pe==0){printf("# time step failed, decreasing time step to: %5.3e\n", flame_params.dt_);}
      }
    }
    if(my_pe==0){printf("# end pseudo time-stepping\n");}
  }

  // Solve system
  flame_params.pseudo_unsteady_ = false;
  flag = KINSetPrintLevel(kinsol_ptr, 1);
  flag = KINSol(kinsol_ptr,
		flame_state,
		KIN_NONE,
		scaler,
		scaler);

  KINGetFuncNorm(kinsol_ptr, &fnorm);

  if((flag==0 || flag==1) && fnorm != 0.0){
    // Get KINSOL stats
    flag = KINGetNumFuncEvals(kinsol_ptr,&nfevals);
    flag = KINGetNumNonlinSolvIters(kinsol_ptr, &nsteps);
#if defined SUNDIALS2 || defined SUNDIALS3
    flag = KINSpilsGetNumPrecEvals(kinsol_ptr,&njacsetups);
    flag = KINSpilsGetNumPrecSolves(kinsol_ptr,&njacsolves);
    flag = KINSpilsGetNumLinIters(kinsol_ptr,&nliniters);
#elif defined SUNDIALS4
    flag = KINGetNumPrecEvals(kinsol_ptr,&njacsetups);
    flag = KINGetNumPrecSolves(kinsol_ptr,&njacsolves);
    flag = KINGetNumLinIters(kinsol_ptr,&nliniters);
#endif
    flag = KINGetFuncNorm(kinsol_ptr, &fnorm);
    flag = KINGetStepLength(kinsol_ptr, &stepnorm);

    if (my_pe==0) {
      printf("Scaled norm of F: %14.7e\n",fnorm);
      printf("Scaled norm of the step: %14.7e\n", stepnorm);
      printf("Number of function evaluations: %ld\n", nfevals);
      printf("Number of nonlinear iterations: %ld\n", nsteps);
      printf("Number of linear iterations: %ld\n", nliniters);
      printf("Number of preconditioner evaluations: %ld\n", njacsetups);
      printf("Number of preconditioner solves: %ld\n", njacsolves);
    }

    // Compute T-Twall
    for(int j=0; j<num_local_points; ++j) {
      int jglobal = j + my_pe*num_local_points;
      temperature_jump[j] =
        (flame_state_ptr[(j+1)*num_states-1]-
         flame_params.wall_temperature_[jglobal])*flame_params.ref_temperature_;

    }
    // Find maximum T-Twall and its location
    max_temperature_jump = FindMaximumParallel(num_local_points,
                                               &flame_params.z_[0],
                                               1,
                                               &temperature_jump[0],
                                               1,
                                               true, // use quadratic
                                               &z_max_temperature_jump);

    // Get min/max of sum(Y_i) and velocity
    min_sum_mass_fraction = minSumMassFractions(flame_state_ptr,flame_params);
    max_sum_mass_fraction = maxSumMassFractions(flame_state_ptr,flame_params);
    min_velocity = minVelocity(flame_state_ptr,flame_params);
    max_velocity = maxVelocity(flame_state_ptr,flame_params);

    // Get maximum of tracked variables
    GetStateMaxima(track_max_state_id,
                   flame_state_ptr,
                   flame_params,
                   &state_maxima,
                   &state_maxima_positions);

    // Print to screen/logfile
    if(my_pe == 0) { //only root prints
      printf("Final  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %6d  %6d  %6d  %6d  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
             current_time,
             flame_params.flame_speed_,
             flame_params.flame_thickness_,
             flame_params.flame_thickness_alpha_,
             z_max_temperature_jump,
             (int)nsteps,
             (int)nfevals,
             (int)njacsetups,
             (int)nliniters,
             0.0,
             0.0,
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

    if(num_prints % flame_params.parser_->field_dt_multiplier() == 0) {
      WriteFieldParallel(1.0,
                         &flame_state_ptr[0],
                         flame_params);
    }
  } // if(flag==0)

  // For error cases
  if(flag<0) {
    if(my_pe == 0) { //only root prints
      printf("Final  0  0  0\n");
    }
  }

  loop_time = getHighResolutionTime() - clock_time;

  // Clear before sens analysis?
  KINFree(&kinsol_ptr);
#if defined SUNDIALS3 || defined SUNDIALS4
  SUNLinSolFree(LS);
#endif
  /*------------------------------------------------------------------------*/
  if(flame_params.flame_speed_sensitivity_) {
    // Perform flame speed sensitivity analysis
    if(my_pe==0) printf("Performing flame speed sensitivity analysis\n");

    clock_time = getHighResolutionTime();
    double multiplier = 2.0; //read from input file?
    rxnSens_t *rxnSensList;
    rxnSensList = new rxnSens_t[num_reactions];

    // 1) Save current/original flame state
    std::vector<double> flame_state_orig;
    flame_state_orig.assign(num_local_states, 0.0); // DUMB
    for(int j=0; j<num_local_states; j++)
      flame_state_orig[j] = flame_state_ptr[j];

    double flame_speed, flame_speed_orig;
    flame_speed_orig = flame_params.flame_speed_;

    // 2) Loop over reactions
    for(int k=0; k<num_reactions; k++) {
      // 2.1) Perturb A factor
      int fwdId = flame_params.mechanism_->getStepIdxOfRxn(k,1);
      int revId = flame_params.mechanism_->getStepIdxOfRxn(k,-1);

      flame_params.reactor_->SetAMultiplierOfStepId(fwdId, multiplier);
      if(revId >= 0 && revId < num_steps)
        flame_params.reactor_->SetAMultiplierOfStepId(revId, multiplier);

      // 2.2) Compute solution
      kinsol_ptr = KINCreate();
      flag = KINInit(kinsol_ptr, ConstPressureFlame, flame_state);
      flag = KINSetErrHandlerFn(kinsol_ptr, ErrorFunction, &flame_params);
      flag = KINSetUserData(kinsol_ptr, &flame_params);
      flag = KINSetFuncNormTol(kinsol_ptr, flame_params.parser_->rel_tol());
      flag = KINSetScaledStepTol(kinsol_ptr, flame_params.parser_->abs_tol());
#if SUNDIALS2
      flag = KINSpgmr(kinsol_ptr, maxl);
      flag = KINSpilsSetMaxRestarts(kinsol_ptr, maxlrst);
      if (flame_params.integrator_type_ == 2) {
        flag = KINSpilsSetPreconditioner(kinsol_ptr, ReactorBBDSetup, ReactorBBDSolve);
      } else if(flame_params.integrator_type_ == 3){
        flag = KINSpilsSetPreconditioner(kinsol_ptr, ReactorAFSetup, ReactorAFSolve);
      } else {
        printf("integrator_type == %d not supported\n",
               flame_params.parser_->integrator_type());
        flag = -1;
      }
#elif SUNDIALS3
      LS = SUNSPGMR(flame_state, PREC_RIGHT, maxl);
      flag = KINSpilsSetLinearSolver(kinsol_ptr, LS);
      flag = SUNSPGMRSetMaxRestarts(LS, maxlrst);
      if (flame_params.integrator_type_ == 2) {
        flag = KINSpilsSetPreconditioner(kinsol_ptr, ReactorBBDSetup, ReactorBBDSolve);
      } else if(flame_params.integrator_type_ == 3){
        flag = KINSpilsSetPreconditioner(kinsol_ptr, ReactorAFSetup, ReactorAFSolve);
      } else {
        printf("integrator_type == %d not supported\n",
               flame_params.parser_->integrator_type());
        flag = -1;
      }
#elif SUNDIALS4
      LS = SUNLinSol_SPGMR(flame_state, PREC_RIGHT, maxl);
      flag = KINSetLinearSolver(kinsol_ptr, LS, NULL);
      flag = SUNLinSol_SPGMRSetMaxRestarts(LS, maxlrst);
      if (flame_params.integrator_type_ == 2) {
	flag = KINSetPreconditioner(kinsol_ptr, ReactorBBDSetup, ReactorBBDSolve);
      } else if(flame_params.integrator_type_ == 3){
	flag = KINSetPreconditioner(kinsol_ptr, ReactorAFSetup, ReactorAFSolve);
      } else {
        printf("integrator_type == %d not supported\n",
               flame_params.parser_->integrator_type());
        flag = -1;
      }
#endif
      flag = KINSetNumMaxIters(kinsol_ptr, maxiter);
      flag = KINSetMaxSetupCalls(kinsol_ptr, mset);
      flag = KINSetPrintLevel(kinsol_ptr, 0);



      flag = KINSol(kinsol_ptr,
                    flame_state,
                    KIN_NONE,
                    scaler,
                    scaler);

      KINFree(&kinsol_ptr);
#if defined SUNDIALS3 || defined SUNDIALS4
      SUNLinSolFree(LS);
#endif

      // 2.3) Get flame speed and compute sensitivity
      flame_speed = flame_params.flame_speed_;
      rxnSensList[k].rxnId = k;
      rxnSensList[k].relSens = fabs(flame_speed - flame_speed_orig)/fabs(flame_speed_orig);
      if(my_pe==0) printf("reaction %d / %d, %s, S: %14.7e\n",
                          k, num_reactions,
	                  flame_params.mechanism_->getReactionName(k),
                          rxnSensList[k].relSens);

      // 2.4) Set A and solution back to original
      flame_params.reactor_->SetAMultiplierOfStepId(fwdId, 1.0);
      if(revId >= 0 && revId < num_steps)
        flame_params.reactor_->SetAMultiplierOfStepId(revId, 1.0);

      for(int j=0; j<num_local_states; j++)
	flame_state_ptr[j] = flame_state_orig[j];

    } // for k<num_reactions

    // Sort sensitivity coefficients in descending order
    qsort((void *)&rxnSensList[0], num_reactions, sizeof(rxnSens_t), compare_rxnSens_t);

    if(my_pe == 0) {
      FILE * sensFile;
      sensFile = fopen("sensAnal_SL","w");
      fprintf(sensFile,"#reac_id  reaction         Srel\n");
      for(int j=0; j<num_reactions; j++) {
        fprintf(sensFile,"%d  %s  %14.7e\n",
                rxnSensList[j].rxnId,
                flame_params.mechanism_->getReactionName(rxnSensList[j].rxnId),
                rxnSensList[j].relSens);
      }
      fclose(sensFile);
    }
    sensanal_time = getHighResolutionTime() - clock_time;


  } //if flame_speed_sensitivity
  /*------------------------------------------------------------------------*/

  if(flame_params.sensitivity_) {
    if(my_pe==0) printf("Computing reaction sensitivity of state variables\n");
      SensitivityAnalysis(flame_state,
                          scaler,
                          &flame_params);
  }

  /*------------------------------------------------------------------------*/

  N_VDestroy_Parallel(flame_state);
  N_VDestroy_Parallel(scaler);

  if(my_pe==0) {
    printf("# Simulation setup time   [s]: %12.5e\n",setup_time);
    printf("# Time in integrator loop [s]: %12.5e\n",loop_time);
    if(flame_params.flame_speed_sensitivity_) printf("# Time in sensitivity loop [s]: %12.5e\n",sensanal_time);
  }
  flame_params.logger_->PrintF(
    "# Simulation setup time   [s]: %12.5e\n",setup_time);
  flame_params.logger_->PrintF(
    "# Time in integrator loop [s]: %12.5e\n",loop_time);

  if(flame_params.integrator_type_ == 2) {
    if(flame_params.superlu_serial_) {
      flame_params.sparse_matrix_->SparseMatrixClean();
    } else {
      flame_params.sparse_matrix_dist_->SparseMatrixClean_dist();
    }
    MPI_Finalize();
  } else {
    MPI_Finalize();
  }
  return 0;
}


// Find maximum and its location
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

// Write state variables to binary file
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
    // Write left BC data
    if(j==num_reactor_states-2) {
      buffer[0] = params.inlet_relative_volume_;
    } else if (j==num_reactor_states-1) {
      buffer[0] = params.inlet_temperature_*params.ref_temperature_;
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
      if(j==num_reactor_states-2){
	buffer[k] = params.rel_vol_[k];
      }
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
  buffer[0] = state[num_reactor_states-2]; //?
  disp = 2*sizeof(int) + sizeof(double) + num_reactor_states*sizeof(char)*64
    + num_reactor_states*(npes*num_local_points+1)*sizeof(double);
  MPI_File_set_view(output_file, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
  MPI_File_write_all(output_file, &buffer[0], 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

  // Interior
  for (int k=0; k<num_local_points; ++k) {
    buffer[k] = state[k*num_reactor_states + num_reactor_states-2];
  }
  disp = 2*sizeof(int) + sizeof(double) + num_reactor_states*sizeof(char)*64
    + sizeof(double)
    + num_reactor_states*(npes*num_local_points+1)*sizeof(double);
  MPI_File_set_view(output_file, disp, MPI_DOUBLE, localarray, "native", MPI_INFO_NULL);
  MPI_File_write_all(output_file, &buffer[0], num_local_points, MPI_DOUBLE, MPI_STATUS_IGNORE);

  MPI_File_close(&output_file);

  MPI_Type_free(&localarray);
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
    velocity = params.rel_vol_[k]*state[(k+1)*num_reactor_states-2];
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
    velocity = params.rel_vol_[k]*state[(k+1)*num_reactor_states-2];
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

    state_max = FindMaximumParallel(num_local_points,
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
 *            NULL pointer   */

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
