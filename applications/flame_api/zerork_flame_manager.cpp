
#include <vector>
#include <string>
#include <map>

#ifdef ZERORK_MPI
#include <mpi.h>
#endif

#include "flame_params.h"
#include "kinsol_functions.h"

#ifdef ZERORK_MPI
#include <nvector/nvector_parallel.h>
#else
#include <nvector/nvector_serial.h>
#endif
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

#include "zerork_flame_manager.h"

#include "zerork/utilities.h"
#include "file_utilities.h"

#include "ZeroRKFlameAPITesterIFP.h"


namespace {

void kin_err_handler(int error_code, const char *module,
                     const char *function, char *msg, void *user_data) {

  FlameParams* flame_params = (FlameParams*) user_data;
  int verbosity = 0;
  flame_params->GetIntOption("verbosity",&verbosity);
  if(verbosity > 0 && flame_params->my_pe_ == 0) {
    fprintf(stderr, "[%s ERROR] %s\n", module, function);
    fprintf(stderr, "  %s\n", msg);
  }
  flame_params->num_kinsol_errors_ += 1;
}

} //anonymous namespace

ZeroRKFlameManager::ZeroRKFlameManager() {
  reactor_   = nullptr;
  transport_ = nullptr;

  tried_init_  = false;
  init_status_ = ZERORK_FLAME_STATUS_SUCCESS;

  //Default options
  int_options_["verbosity"] = 0;
  double_options_["reference_temperature"] = 4000.0;
  string_options_["mechanism_filename"] = std::string("mech.dat");
  string_options_["therm_filename"] = std::string("therm.dat");
  string_options_["transport_filename"] = std::string("tran.dat");
  string_options_["transport_model"] = "ConstantLewis";
  string_options_["constant_lewis_setting"] = "Tfix"; //GridPoint, Tfix, Default
  string_options_["temperature_fix_setting"] = "InOutMix"; //Absolute, InOutMix
  string_options_["kinsol_strategy"] = "NONE"; //"LINESEARCH"
  int_options_["constant_lewis_grid_point"] = 0;
  double_options_["temperature_fix_value"] = 0.85;

  //Solver options
  int_options_["integrator_type"] = 3;
  int_options_["store_jacobian"]  = 1;
  int_options_["convective_scheme_type"] = 0;
  int_options_["max_subiterations"] = 10;
  int_options_["steady_max_iterations"] = 100;
  double_options_["relative_tolerance"] = 1.0e-3;
  double_options_["absolute_tolerance"] = 1.0e-20;
  double_options_["step_limiter"] = 1.0e300;
  int_options_["pseudo_unsteady"] = 1;
  double_options_["pseudo_unsteady_dt"] = 1.0e-7;
  double_options_["pseudo_unsteady_min_dt"] = 1.0e-10;
  double_options_["pseudo_unsteady_max_dt"] = 5.0e-2;
  int_options_["pseudo_unsteady_max_iterations"] = 20;
  double_options_["pseudo_unsteady_time"] = 1.0;

  //File-output Options
  string_options_["mechanism_parsing_log_filename"] = std::string(zerork::utilities::null_filename);
  string_options_["transport_parsing_log_filename"] = std::string(zerork::utilities::null_filename);
}

zerork_flame_status_t ZeroRKFlameManager::ReadOptionsFile(const std::string& options_filename) {
  std::unique_ptr<ZeroRKFlameAPITesterIFP> inputFileDBptr;
  try {
    inputFileDBptr = std::make_unique<ZeroRKFlameAPITesterIFP>(options_filename);
  } catch (const std::runtime_error& e) {
    return ZERORK_FLAME_STATUS_FAILED_OPTIONS_PARSE;
  }
  const ZeroRKFlameAPITesterIFP& inputFileDB(*inputFileDBptr);

  //TODO: Don't over-ride if options file used default value?
  int_options_["verbosity"] = inputFileDB.verbosity();
  int_options_["steady_max_iterations"] = inputFileDB.steady_max_iterations();
  int_options_["max_subiterations"] = inputFileDB.max_subiterations();
  int_options_["integrator_type"] = inputFileDB.integrator_type();
  int_options_["store_jacobian"]  = inputFileDB.store_jacobian();
  int_options_["convective_scheme_type"] = inputFileDB.convective_scheme_type();
  int_options_["pseudo_unsteady"] = inputFileDB.pseudo_unsteady();
  int_options_["pseudo_unsteady_max_iterations"] = inputFileDB.pseudo_unsteady_max_iterations();
  int_options_["constant_lewis_grid_point"] = inputFileDB.constant_lewis_grid_point();

  double_options_["absolute_tolerance"] = inputFileDB.absolute_tolerance();
  double_options_["relative_tolerance"] = inputFileDB.relative_tolerance();
  double_options_["reference_temperature"] = inputFileDB.reference_temperature();
  double_options_["pseudo_unsteady_dt"] = inputFileDB.pseudo_unsteady_dt();
  double_options_["pseudo_unsteady_min_dt"] = inputFileDB.pseudo_unsteady_min_dt();
  double_options_["pseudo_unsteady_max_dt"] = inputFileDB.pseudo_unsteady_max_dt();
  double_options_["pseudo_unsteady_time"] = inputFileDB.pseudo_unsteady_time();
  double_options_["step_limiter"] = inputFileDB.step_limiter();
  double_options_["temperature_fix_value"] = inputFileDB.temperature_fix_value();

  string_options_["mechanism_parsing_log_filename"] = inputFileDB.mechanism_parsing_log();
  string_options_["transport_parsing_log_filename"] = inputFileDB.transport_parsing_log();
  string_options_["mechanism_filename"] = inputFileDB.mechanism_file();
  string_options_["therm_filename"] = inputFileDB.therm_file();
  string_options_["transport_filename"] = inputFileDB.transport_file();
  string_options_["transport_model"] = inputFileDB.transport_model();
  string_options_["constant_lewis_setting"] = inputFileDB.constant_lewis_setting();
  string_options_["temperature_fix_setting"] = inputFileDB.temperature_fix_setting();
  string_options_["kinsol_strategy"] = inputFileDB.kinsol_strategy();
  return ZERORK_FLAME_STATUS_SUCCESS;
}

zerork_flame_status_t ZeroRKFlameManager::LoadMechanism() {
  try {
    reactor_ = std::make_shared<ConstPressureReactor>(string_options_["mechanism_filename"].c_str(),
                                        string_options_["therm_filename"].c_str(),
                                        string_options_["mechanism_parsing_log_filename"].c_str(),
                                        COMPRESSED_COL_STORAGE,
                                        101325.0);
  } catch (const std::runtime_error& e) {
    reactor_ = nullptr;
    return ZERORK_FLAME_STATUS_FAILED_MECHANISM_PARSE;
  }

  // setup transport interface
  transport_ = std::shared_ptr<transport::MassTransportInterface>(
                  transport::InterfaceFactory::CreateMassBased(string_options_["transport_model"])
               );

  std::vector<std::string> transport_files;
  transport_files.push_back(string_options_["transport_filename"]);

  int error_code = transport_->Initialize(reactor_->GetMechanism(),
                                  transport_files,
                                  string_options_["transport_parsing_log_filename"]);
  if(error_code != transport::NO_ERROR) {
    printf("# ERROR: Could not Initialize MassTransportInterface for files:\n"
           "#            mechanism      file = %s\n"
           "#            thermodynamics file = %s\n"
           "#            transport      file = %s\n"
           "#            log            file = %s\n"
           "#        Initialize returned error code = %d\n",
           string_options_["mechanism_filename"].c_str(),
           string_options_["therm_filename"].c_str(),
           string_options_["transport_filename"].c_str(),
           string_options_["transport_parsing_log_filename"].c_str(),
           error_code);
    return ZERORK_FLAME_STATUS_FAILED_TRANSPORT_PARSE;
  }
  return ZERORK_FLAME_STATUS_SUCCESS;
}

zerork_flame_status_t ZeroRKFlameManager::FinishInit() {
  if(!tried_init_) {
    tried_init_ = true;
 
    //parse mechanism if not yet
    if(reactor_ == nullptr) {
      init_status_ = this->LoadMechanism();
    }
  }
  return init_status_;
}

zerork_flame_status_t ZeroRKFlameManager::Solve(int num_grid_points, const double* grid_points, double P,
                                                double* flame_speed, double* T, double* mass_fractions)
{
  const int num_species     = reactor_->GetNumSpecies();
  const int integrator_type = int_options_["integrator_type"];

  std::vector<double> grid(grid_points, grid_points+num_grid_points);

  reactor_->SetReferenceTemperature(double_options_["reference_temperature"]);
  reactor_->SetPressure(P);

  // Initialize flame params
#ifdef ZERORK_MPI
  FlameParams flame_params(reactor_.get(), transport_.get(), grid, *flame_speed, T, mass_fractions, P, MPI_COMM_WORLD, *this);
#else
  FlameParams flame_params(reactor_.get(), transport_.get(), grid, *flame_speed, T, mass_fractions, P, *this);
#endif
  if(flame_params.error_status_ != 0) {
    // We may want to add some status'es for more information on what/why FlameParams errored
    return ZERORK_FLAME_STATUS_FAILED_SOLVE;
  }

  //Declare variables
  int maxl, maxlrst;

  long int num_local_states = flame_params.num_local_points_*reactor_->GetNumStates();
  long int total_states = flame_params.z_.size()*reactor_->GetNumStates();

  // KINSol integrator Stats
  long int nsteps, nfevals, nliniters, njacsetups, njacsolves;

  // Initialize state variables/scalers vectors and pointers
#ifdef ZERORK_MPI
  N_Vector flame_state = N_VMake_Parallel(flame_params.comm_, num_local_states, total_states, flame_params.y_.data());
  N_Vector flame_state_old = N_VMake_Parallel(flame_params.comm_, num_local_states, total_states, flame_params.y_old_.data());
  N_Vector scaler      = N_VNew_Parallel(flame_params.comm_, num_local_states, total_states);
#else
  N_Vector flame_state = N_VMake_Serial(total_states, flame_params.y_.data());
  N_Vector flame_state_old = N_VMake_Serial(total_states, flame_params.y_old_.data());
  N_Vector scaler      = N_VNew_Serial(total_states);
#endif

  // Initialize scaler
  N_VConst(1.0, scaler);

  //----------------------------------------------------------------------------
  // Setup KINSOL solver
  // Create KINSOL pointer
  // KINSOL memory pointer and linear solver
  void* kinsol_ptr = KINCreate();

  // Initialize KINSOL module with RHS function and state vector
  int kinsol_flag = KINInit(kinsol_ptr, ConstPressureFlame, flame_state);

  // Set user data
  kinsol_flag = KINSetUserData(kinsol_ptr, &flame_params);
  kinsol_flag = KINSetErrHandlerFn(kinsol_ptr, &kin_err_handler, &flame_params);

  // Set tolerances
  // RHS(y) < fnormtol
  kinsol_flag = KINSetFuncNormTol(kinsol_ptr, double_options_["relative_tolerance"]);
  // Step tolerance: dy < steptol
  kinsol_flag = KINSetScaledStepTol(kinsol_ptr, double_options_["absolute_tolerance"]);

  // Setup KINSol
  maxl = 1000; // max Krylov dimension
  maxlrst = 0; // max linear solver restarts

#ifdef SUNDIALS2
  // Initialize Linear Solver
  kinsol_flag = KINSpgmr(kinsol_ptr, maxl);
  kinsol_flag = KINSpilsSetMaxRestarts(kinsol_ptr, maxlrst);
  // Set preconditioner
  if (integrator_type == 2) {
    kinsol_flag = KINSpilsSetPreconditioner(kinsol_ptr, ReactorBBDSetup, ReactorBBDSolve);
  } else if(integrator_type == 3){
    kinsol_flag = KINSpilsSetPreconditioner(kinsol_ptr, ReactorAFSetup, ReactorAFSolve);
  }
#elif SUNDIALS3
  // Initialize Linear Solver
  SUNLinearSolver LS = SUNSPGMR(flame_state, PREC_RIGHT, maxl);
  kinsol_flag = KINSpilsSetLinearSolver(kinsol_ptr, LS);
  kinsol_flag = SUNSPGMRSetMaxRestarts(LS, maxlrst);

  // Set preconditioner
  if (integrator_type == 2) {
    kinsol_flag = KINSpilsSetPreconditioner(kinsol_ptr, ReactorBBDSetup, ReactorBBDSolve);
  } else if(integrator_type == 3){
    kinsol_flag = KINSpilsSetPreconditioner(kinsol_ptr, ReactorAFSetup, ReactorAFSolve);
  }
#elif SUNDIALS4
  // Initialize Linear Solver
  SUNLinearSolver LS = SUNLinSol_SPGMR(flame_state, PREC_RIGHT, maxl);
  kinsol_flag = KINSetLinearSolver(kinsol_ptr, LS, NULL);
  kinsol_flag = SUNLinSol_SPGMRSetMaxRestarts(LS, maxlrst);

  // Set preconditioner
  if (integrator_type == 2) {
    kinsol_flag = KINSetPreconditioner(kinsol_ptr, ReactorBBDSetup, ReactorBBDSolve);
  } else if(integrator_type == 3){
    kinsol_flag = KINSetPreconditioner(kinsol_ptr, ReactorAFSetup, ReactorAFSolve);
  }
#endif

  //0 for default, 1 for exact Newton, > 1 for modified Newton
  kinsol_flag = KINSetMaxSetupCalls(kinsol_ptr, int_options_["max_subiterations"]);

  if(flame_params.my_pe_ == 0) {
    // 0 for no info, 1 for scaled l2 norm, 3 for additional linear solver info
    int print_level = std::max(0, int_options_["verbosity"] - 1);
    kinsol_flag = KINSetPrintLevel(kinsol_ptr, print_level);
  } else {
    kinsol_flag = KINSetPrintLevel(kinsol_ptr, 0);
  }

  //----------------------------------------------------------------------------

  // Get time
  double clock_time = zerork::getHighResolutionTime();
  double fnorm, stepnorm;

  int strategy = KIN_NONE;
  if( string_options_["kinsol_strategy"] == "LINESEARCH" ) {
    strategy = KIN_LINESEARCH;
  }

  // Try to solve system with pure steady RHS
  flame_params.pseudo_unsteady_ = false;

  kinsol_flag = KINSetNumMaxIters(kinsol_ptr, int_options_["steady_max_iterations"]);

  N_VScale(1.0, flame_state, flame_state_old);

  if(int_options_["verbosity"] > 0 && flame_params.my_pe_ == 0) {
    printf("Starting steady solve\n");
  }
  kinsol_flag = KINSol(kinsol_ptr,
                       flame_state,
                       strategy,
                       scaler,
                       scaler);

  // Kinsol error
  if(kinsol_flag<0) {
    if(int_options_["verbosity"] > 0 && flame_params.my_pe_ == 0) {
      printf("Failed steady solve\n");
    }

    strategy = KIN_LINESEARCH;

    N_VScale(1.0, flame_state_old, flame_state);
    if(flame_params.my_pe_ == 0) {
      flame_params.flame_speed_ = flame_params.y_[num_species]*flame_params.inlet_relative_volume_;
    }
#ifdef ZERORK_MPI
    MPI_Bcast(&flame_params.flame_speed_, 1, MPI_DOUBLE, 0, flame_params.comm_);
#endif

    flame_params.pseudo_unsteady_ = int_options_["pseudo_unsteady"] != 0;
    // Pseudo-unsteady
    if(flame_params.pseudo_unsteady_) {
      kinsol_flag = KINSetNumMaxIters(kinsol_ptr, int_options_["pseudo_unsteady_max_iterations"]);
      double pseudo_time = 0.0;
      int pseudo_successful_steps = 0;
      double kinstart_time, kinend_time;
      flame_params.dt_ = double_options_["pseudo_unsteady_dt"];

      while(pseudo_time < double_options_["pseudo_unsteady_time"]) {
        N_VScale(1.0, flame_state, flame_state_old);
        double flame_speed_old = flame_params.flame_speed_;

        // Solve system
        kinstart_time = zerork::getHighResolutionTime();
        if(int_options_["verbosity"] > 0 && flame_params.my_pe_ == 0) {
          printf("Starting pseudo_unsteady solve at time %g s with step %g s (flame speed %g)\n", pseudo_time, flame_params.dt_, flame_params.flame_speed_);
        }
        kinsol_flag = KINSol(kinsol_ptr,
                      flame_state,
                      strategy,
                      scaler,
                      scaler);
        kinend_time = zerork::getHighResolutionTime();

        KINGetNumFuncEvals(kinsol_ptr,&nfevals);
        KINGetFuncNorm(kinsol_ptr, &fnorm);

        if(pseudo_successful_steps>  5 &&
           fabs((flame_speed_old-flame_params.flame_speed_)/flame_params.flame_speed_) < 1.0e-5 &&
           flame_params.dt_ > 1.0e-6) {
          break;
        }
        if((kinsol_flag==0 || kinsol_flag==1) && fnorm != 0.0){
          // Success
          pseudo_time += flame_params.dt_;
          pseudo_successful_steps += 1;

          // Increase time step if it's converging quickly
          if(nfevals <= 5) {
            flame_params.dt_ *= 2.0;
          } else if (nfevals <= 10) {
            flame_params.dt_ *= 1.2;
          }
          flame_params.dt_ = std::min(flame_params.dt_,double_options_["pseudo_unsteady_max_dt"]);
        } else {
          // Failure, reset state and decrease timestep
          if(int_options_["verbosity"] > 0 && flame_params.my_pe_ == 0) {
            printf("Failed pseudo_unsteady solve at time %g s (flame speed %g)\n", pseudo_time, flame_params.flame_speed_);
          }
          N_VScale(1.0, flame_state_old, flame_state);
          flame_params.flame_speed_ = flame_speed_old;
          flame_params.dt_ *= 0.5;
          pseudo_successful_steps = 0;
          if(flame_params.dt_ < double_options_["pseudo_unsteady_min_dt"]){
            return ZERORK_FLAME_STATUS_FAILED_SOLVE;
          }
        }
      }
    }

    // Solve system with pure steady RHS
    flame_params.pseudo_unsteady_ = false;

    kinsol_flag = KINSetNumMaxIters(kinsol_ptr, int_options_["steady_max_iterations"]);

    if(int_options_["verbosity"] > 0 && flame_params.my_pe_ == 0) {
      printf("Starting steady solve\n");
    }
    kinsol_flag = KINSol(kinsol_ptr,
                         flame_state,
                         strategy,
                         scaler,
                         scaler);

    if(kinsol_flag < 0 && strategy != KIN_NONE) {
      //sometimes linesearch fails when KIN_NONE succeeds
      kinsol_flag = KINSol(kinsol_ptr,
                           flame_state,
                           KIN_NONE,
                           scaler,
                           scaler);
    }
  }

  KINFree(&kinsol_ptr);
#if defined SUNDIALS3 || defined SUNDIALS4
  SUNLinSolFree(LS);
#endif

  N_VDestroy(flame_state);
  N_VDestroy(scaler);
  N_VDestroy(flame_state_old);

  // Kinsol error
  if(kinsol_flag<0) {
    if(int_options_["verbosity"] > 0 && flame_params.my_pe_ == 0) {
      printf("Failed steady solve\n");
    }
    return ZERORK_FLAME_STATUS_FAILED_SOLVE;
  } else {

    *flame_speed = flame_params.flame_speed_;
    flame_params.GetTemperatureAndMassFractions(T, mass_fractions);

    if(int_options_["verbosity"] > 0 && flame_params.my_pe_ == 0) {
      printf("Succeeded steady solve\n");
    }
    return ZERORK_FLAME_STATUS_SUCCESS;
  }
}

