#include <stdlib.h>
#include <stdio.h>

#include <string>
#include <vector>

#include <cvodes/cvodes.h>            // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros

#ifdef SUNDIALS2
#include <cvode/cvode_dense.h>      // prototypes for CVDense
#include <cvode/cvode_spgmr.h>      // prototypes & constants for CVSPGMR
#elif SUNDIALS3
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <cvode/cvode_spils.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#elif SUNDIALS4
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif

#include "integrated_function.h"
#include "user_functions_gsa.h"
#include "complete_solver_gsa.h"

using zerork::getHighResolutionTime;

const int MAX_LINE_LENGTH = 8192;

// FindMaximum assumes x is in increasing order
static double FindMaximum(const size_t num_points,
                          const double x[],
                          const size_t x_stride,
                          const double f[],
                          const size_t f_stride,
                          const bool use_quadratic,
                          double *x_at_max);

static double QuadraticMaximum(const double x[],
                               const size_t x_stride,
                               const double f[],
                               const size_t f_stride,
                               double *x_at_max);

void SetupCompleteSolver(void *cvode_memory,
                     aux_cvode_structs_t &aux_cvode,
                         UserData *user_data,
                         N_Vector cvode_state)
{
  int error_flag;
  //user_data->GetInitialState(NV_DATA_S(cvode_state));

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  error_flag = CVodeInit(cvode_memory,
                         VariableVolumeRHS,
                         user_data->GetInitialTime(),
                         cvode_state);
  if (CheckFlag(&error_flag, "CVodeInit", 1)) exit(-1);


  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  error_flag = CVodeSStolerances(cvode_memory,
                                 user_data->GetParser()->relTol(),
                                 user_data->GetParser()->absTol());
  if (CheckFlag(&error_flag, "CVodeSStolerances", 1)) exit(-1);

  /* Set the pointer to user-defined data */
  error_flag = CVodeSetUserData(cvode_memory,
                                user_data);
  if(CheckFlag(&error_flag, "CVodeSetUserData", 1)) exit(-1);

  /* Set the quadrature function for the chemical heat release rate */
  error_flag = CVodeQuadInit(cvode_memory,
                             ChemicalHRRIntegrand,
                             user_data->quadrature_);
  if(CheckFlag(&error_flag, "CVodeQuadInit", 1)) exit(-1);
  error_flag = CVodeQuadSStolerances(cvode_memory,
                                     user_data->quadrature_rtol_,
                                     user_data->quadrature_atol_);
  if(CheckFlag(&error_flag, "CVodeQuadSStolerances", 1)) exit(-1);

  error_flag = CVodeSetQuadErrCon(cvode_memory,
                                  user_data->quadrature_controls_step_);
  if(CheckFlag(&error_flag, "CVodeSetQuadErrCons", 1)) exit(-1);


  /* Set the maximum number of internal steps per CVode call and the maximum
   * allowable internal steps. */
  error_flag = CVodeSetMaxOrd(cvode_memory,
                              user_data->GetParser()->maxOrder());
  if (CheckFlag(&error_flag, "CVodeSetMaxOrd", 1)) exit(-1);
  error_flag = CVodeSetMaxNumSteps(cvode_memory,
                                   user_data->GetParser()->maxSteps());
  if (CheckFlag(&error_flag, "CVodeSetMaxNumSteps", 1)) exit(-1);

  error_flag = CVodeSetMaxStep(cvode_memory,
                               user_data->GetParser()->maxDtInternal());
  if (CheckFlag(&error_flag, "CVodeSetMaxStep", 1)) exit(-1);

  error_flag = CVodeSetNonlinConvCoef(cvode_memory,
                                      user_data->GetParser()->cvNlConvCoeff());
  if (CheckFlag(&error_flag, "CVodeSetNonlinConvCoef", 1)) exit(-1);

  // Attach the linear solver
  //SetupDenseSolver(cvode_memory,
  //                 user_data);
  SetupSparseSolver(cvode_memory,
                      aux_cvode,
                    user_data,
                    cvode_state);
}

void SetupSparseSolver(void *cvode_memory,
                       aux_cvode_structs_t &aux_cvode,
                       UserData *user_data,
                       N_Vector cvode_state)
{
  int flag;
#ifdef SUNDIALS2
  /* Call CVSpgmr to specify the linear solver CVSPGMR
     with left preconditioning and the maximum Krylov dimension maxl */
  flag = CVSpgmr(cvode_memory,
                 PREC_LEFT,
                 user_data->GetParser()->maxKrylovDim());
  if(CheckFlag(&flag, "CVSpgmr", 1)) exit(-1);

  flag = CVSpilsSetGSType(cvode_memory, MODIFIED_GS);
  if(CheckFlag(&flag, "CVSpilsSetGSType", 1)) exit(-1);

  /* Set preconditioner setup and solve routines Precond and PSolve,
     and the pointer to the user-defined block data */
  flag = CVSpilsSetPreconditioner(cvode_memory,
                                  VariableVolumePreconditionerSetup,
  				  VariableVolumePreconditionerSolve);
  if(CheckFlag(&flag, "CVSpilsSetPreconditioner", 1)) exit(-1);
  flag = CVSpilsSetEpsLin(cvode_memory,
                          user_data->GetParser()->cvEpsLin());
  if (CheckFlag(&flag, "CVSpilsSetEpsLin", 1)) exit(-1);
#elif SUNDIALS3
  aux_cvode.alloc = 2;
  aux_cvode.LS = SUNSPGMR(cvode_state, PREC_LEFT, 5); //!?!?!
  flag = CVSpilsSetLinearSolver(cvode_memory, aux_cvode.LS);
  if(CheckFlag(&flag, "CVSpilsSetLinearSolver", 1)) exit(-1);
  flag = SUNSPGMRSetGSType(aux_cvode.LS, MODIFIED_GS);
  if(CheckFlag(&flag, "SUNSPGMRSetGSType", 1)) exit(-1);
  flag = CVSpilsSetPreconditioner(cvode_memory,
                                  VariableVolumePreconditionerSetup,
                                  VariableVolumePreconditionerSolve);
  if(CheckFlag(&flag, "CVSpilsSetPreconditioner", 1)) exit(-1);
  flag = CVSpilsSetEpsLin(cvode_memory,
                          user_data->GetParser()->cvEpsLin());
  if (CheckFlag(&flag, "CVSpilsSetEpsLin", 1)) exit(-1);
#elif SUNDIALS4
  aux_cvode.alloc = 2;
  aux_cvode.NLS = SUNNonlinSol_Newton(cvode_state);
  flag = CVodeSetNonlinearSolver(cvode_memory, aux_cvode.NLS);
  if(CheckFlag(&flag, "CVodeSetNonlinearSolver", 1)) exit(-1);
  aux_cvode.LS = SUNLinSol_SPGMR(cvode_state, PREC_LEFT, 5);
  flag = CVodeSetLinearSolver(cvode_memory, aux_cvode.LS, NULL);
  if(CheckFlag(&flag, "CVodeSetLinearSolver", 1)) exit(-1);

  flag = CVodeSetPreconditioner(cvode_memory,
                                VariableVolumePreconditionerSetup,
                                VariableVolumePreconditionerSolve);
  if(CheckFlag(&flag, "CVodeSetPreconditioner", 1)) exit(-1);

  flag = CVodeSetEpsLin(cvode_memory,
                        user_data->GetParser()->cvEpsLin());
  if (CheckFlag(&flag, "CVodeSetEpsLin", 1)) exit(-1);
#endif

  // flag = CVSpilsSetJacTimesVecFn(cvode_mem, sparse_jac_v);
  // if(check_flag(&flag, "CVSpilsSetJacTimesVecFn", 1)) exit(-1);

}
void SetupDenseSolver(void *cvode_memory,
                      aux_cvode_structs_t &aux_cvode,
                      UserData *user_data,
                      N_Vector cvode_state)
{
  int error_flag;
  // Dense solver
#ifdef SUNDIALS2
  error_flag = CVDense(cvode_memory,
                       user_data->GetNumStates());
  if (CheckFlag(&error_flag, "CVDense", 1)) exit(-1);

  // Attach the dense jacobian function
  error_flag = CVDlsSetDenseJacFn(cvode_memory,
                                  VariableVolumeDenseJacobian);
  if (CheckFlag(&error_flag, "CVDlsSetDenseJacFn", 1)) exit(-1);
#elif SUNDIALS3
  aux_cvode.alloc = 1;
  aux_cvode.A = SUNDenseMatrix(user_data->GetNumStates(), user_data->GetNumStates());
  if(CheckFlag((void *)aux_cvode.A, "SUNDenseMatrix", 0)) exit(-1);
  aux_cvode.LS = SUNDenseLinearSolver(cvode_state, aux_cvode.A);
  if(CheckFlag((void *)aux_cvode.LS, "SUNDenseLinearSolver", 0)) exit(-1);
  error_flag = CVDlsSetLinearSolver(cvode_memory, aux_cvode.LS, aux_cvode.A);
  if(CheckFlag(&error_flag, "CVDlsSetLinearSolver", 1)) exit(-1);
  error_flag = CVDlsSetJacFn(cvode_memory, VariableVolumeDenseJacobian);
  if(CheckFlag(&error_flag, "CVDlsSetJacFn", 1)) exit(-1);
#elif SUNDIALS4
  aux_cvode.alloc = 1;
  aux_cvode.A = SUNDenseMatrix(user_data->GetNumStates(), user_data->GetNumStates());
  if(CheckFlag((void *)aux_cvode.A, "SUNDenseMatrix", 0)) exit(-1);
  aux_cvode.LS = SUNLinSol_Dense(cvode_state, aux_cvode.A);
  if(CheckFlag((void *)aux_cvode.LS, "SUNLinSol_Dense", 0)) exit(-1);
  error_flag = CVodeSetLinearSolver(cvode_memory, aux_cvode.LS, aux_cvode.A);
  if(CheckFlag(&error_flag, "CVodeSetLinearSolver", 1)) exit(-1);
  error_flag = CVodeSetJacFn(cvode_memory, VariableVolumeDenseJacobian);
  if(CheckFlag(&error_flag, "CVodeSetJacFn", 1)) exit(-1);
#endif
}


//  const int NUM_FIXED_RESULTS = 6; // from complete_solver_gsa.h
//
//  results[0] = current_time;
//  results[1] = max_dp_dt;
//  results[2] = time_max_dp_dt;
//  results[3] = total_heat_release;
//  results[4] = max heat release rate (piece-wise linear)
//  results[5] = max heat release rate (piece-wise backward quadratic)
//  for(int j=0; j<num_burn_fractions; ++j) {
//    results[j+NUM_FIXED_RESULTS] = burn_fraction_time[j];
//  }
//
//  where
//
//   num_burn_fractions =
//     static_cast<int>(user_data->GetParser()->burnFraction().size());
void DestroyAuxCVodeStructs(aux_cvode_structs_t &aux_cvode)
{
#if defined SUNDIALS3 || defined SUNDIALS4
  if(aux_cvode.alloc > 0) {
    SUNLinSolFree(aux_cvode.LS);
  }
  if(aux_cvode.alloc == 1) {
    SUNMatDestroy(aux_cvode.A);
  }
#endif
#if defined SUNDIALS4
  if(aux_cvode.alloc == 2) {
    SUNNonlinSolFree(aux_cvode.NLS);
  }
#endif
  aux_cvode.alloc = 0;
}


int SolveVariableVolume(void *cvode_memory,
                        UserData *user_data,
                        N_Vector cvode_state,
                        const double a_multipliers[],
                        const TimeHistoryParams &thist_params,
                        std::string *thist,
                        double *results,
                        double *simulation_time)
{
  const int num_reactions = user_data->GetReactor()->GetNumReactions();
  const int num_burn_fractions = static_cast<int>(user_data->GetParser()->burnFraction().size());
  const int num_hrr_powers = static_cast<int>(user_data->num_quadratures_ - 1);
  const int num_results = NUM_FIXED_RESULTS + num_hrr_powers + num_burn_fractions;

  double sim_start = getHighResolutionTime();

  double current_time;
  double final_time;
  int num_steps=0;
  int error_flag;
  double current_pressure;
  double previous_time, previous_pressure;
  double print_time;
  double current_dp_dt, max_dp_dt, time_max_dp_dt;
  std::string print_line;
  double total_heat_release;
  double current_heat_release_rate;
  double inv_reference_heat_release;
  double current_burn_fraction;
  double previous_burn_fraction;
  bool *exceeded_burn_fraction;
  double *burn_fraction_time;
  std::vector<double> time_history;
  std::vector<double> heat_release_rate_history;
  double time_at_history_max;
  double ode_time_period;
  double quadrature_time;

  if(num_burn_fractions > 0) {
    burn_fraction_time = new double[num_burn_fractions];
    exceeded_burn_fraction = new bool[num_burn_fractions];

    for(int j=0; j<num_burn_fractions; ++j) {
      exceeded_burn_fraction[j] = false;
      burn_fraction_time[j]     = 1.0e300;
    }
    inv_reference_heat_release = 1.0/user_data->GetReferenceHeatRelease();
  }
  *simulation_time = 0.0;
  for(int j=0; j<num_results; ++j) {
    results[j] = 0.0;
  }

  // ------------------------------------------------------------------------
  // set the initial time and state, and the final time
  //user_data.GetInitialState(NV_DATA_S(cvode_state));
  current_time = user_data->GetInitialTime();
  final_time = user_data->GetParser()->finalTime();

  // write the header lines if necessary
  if(thist_params.record_history || thist_params.echo_stdout) {
    BuildTimeHistoryHeader(cvode_memory,
                           current_time,
                           cvode_state,
                           user_data,
                           print_line);
    if(thist_params.echo_stdout) {
      printf("%s",print_line.c_str());
    }
    if(thist_params.record_history) {
      thist->clear();
      thist->append(print_line);
    }
  }
  // Set new A-Factor multipliers
  if(a_multipliers != NULL) {
    for(int j=0; j<num_reactions; ++j) {
      user_data->GetReactor()->SetAMultiplierOfForwardReactionId(j,
							     a_multipliers[j]);
      user_data->GetReactor()->SetAMultiplierOfReverseReactionId(j,
							     a_multipliers[j]);
      // nothing is set if reverse reaction doesn't exist
      // SetAMultiplierOf_______ReactionId(...) returns INDEX_OUT_OF_RANGE
      // ReactorError code in this case
    }
  }

  // reset initial state
  user_data->GetInitialState(NV_DATA_S(cvode_state));
  // reset cvode
  error_flag = CVodeReInit(cvode_memory, current_time, cvode_state);
  if(CheckFlag(&error_flag, "CVodeReInit", 1)) {
    return error_flag;
  }
  // reset the quadrature variables
  for(int j=0; j<user_data->num_quadratures_; ++j) {
    NV_Ith_S(user_data->quadrature_,j) = 0.0;
  }
  // reset cvode quadrature
  error_flag = CVodeQuadReInit(cvode_memory, user_data->quadrature_);
  if(CheckFlag(&error_flag, "CVodeQuadReInit", 1)) {
    return error_flag;
  }

  // record initial point
  total_heat_release = 0.0;
  current_burn_fraction = 0.0;
  max_dp_dt      = 0.0;
  time_max_dp_dt = current_time;
  current_pressure =
    user_data->GetPressure(current_time,NV_DATA_S(cvode_state));
  current_heat_release_rate = ChemicalHeatReleaseRate(current_time,
                                                      NV_DATA_S(cvode_state),
                                                      user_data);
  time_history.push_back(current_time);
  heat_release_rate_history.push_back(current_heat_release_rate);

  // set the last print more than one min_time_step in the past
  print_time = current_time - 1.1*thist_params.min_time_step;


  while(current_time < final_time) {
    previous_time = current_time;
    previous_pressure = current_pressure;
    previous_burn_fraction = current_burn_fraction;
    error_flag =  CVode(cvode_memory,
                        final_time,
                        cvode_state,
                        &current_time,
                        CV_ONE_STEP);
    ++num_steps;

    // write the header lines if necessary
    if(thist_params.record_history || thist_params.echo_stdout) {

      if((current_time > thist_params.min_time_step + print_time) ||
         (num_steps%thist_params.step_print_period == 0)) {

        BuildTimeHistoryLine(cvode_memory,
                             current_time,
                             cvode_state,
                             user_data,
                             print_line);

        if(thist_params.echo_stdout) {
          printf("%s",print_line.c_str());
        }
        if(thist_params.record_history) {
          thist->append(print_line);
        }
        print_time = current_time;
      }
    }

    if(CheckFlag(&error_flag, "CVODE ERROR", 1)) {
      // TO DO: add logic to restart the integrator a user-specified number
      //        of times

      // make sure the time history line is recorded in the case
      // of a CVODE failure
      if((thist_params.record_history || thist_params.echo_stdout) &&
         print_time < current_time) {

        BuildTimeHistoryLine(cvode_memory,
                             current_time,
                             cvode_state,
                             user_data,
                             print_line);

        if(thist_params.echo_stdout) {
          printf("%s",print_line.c_str());
        }
        if(thist_params.record_history) {
          thist->append(print_line);
        }
      }
      // output the raw state values and derivatives at the failure
      PrintStateDerivatives(current_time,
                            cvode_state,
                            user_data);

      return error_flag;
    }

    current_pressure =
      user_data->GetPressure(current_time,NV_DATA_S(cvode_state));
    current_dp_dt = (current_pressure-previous_pressure)/
      (current_time-previous_time);
    if(current_dp_dt > max_dp_dt) {
      max_dp_dt = current_dp_dt;
      time_max_dp_dt = current_time;
    }
    // record heat release rate history
    current_heat_release_rate = ChemicalHeatReleaseRate(current_time,
                                                        NV_DATA_S(cvode_state),
                                                        user_data);
    time_history.push_back(current_time);
    heat_release_rate_history.push_back(current_heat_release_rate);
    error_flag = CVodeGetQuad(cvode_memory,
                              &quadrature_time,
                              user_data->quadrature_);
    if(CheckFlag(&error_flag, "CVODE ERROR", 1)) {
      return error_flag;
    }
    total_heat_release = NV_Ith_S(user_data->quadrature_,0);

    if(num_burn_fractions > 0) {
      current_burn_fraction = total_heat_release*inv_reference_heat_release;
      //printf("debug %14.7e  %14.7e  %14.7e\n",
      //       current_time,total_heat_release, current_burn_fraction);
      // scan the burn fraction list to see if any are exceeded on this
      // time step
      for(int j=0; j<num_burn_fractions; ++j) {
        if((current_burn_fraction >=
            user_data->GetParser()->burnFraction()[j]) &&
           exceeded_burn_fraction[j] == false) {
          // record the first time the burn fraction is exceeded

          // TODO: figure out why the first order linear approximation
          //       produces NaNs
          burn_fraction_time[j] =
            current_time - (current_time-previous_time)*
	    (current_burn_fraction-user_data->GetParser()->burnFraction()[j])/
            (current_burn_fraction-previous_burn_fraction);
	  //burn_fraction_time[j] = current_time;
	  exceeded_burn_fraction[j] = true;
        }
      }
    }
  }

  error_flag = CVodeGetQuad(cvode_memory,
                            &quadrature_time,
                            user_data->quadrature_);
  if(CheckFlag(&error_flag, "CVODE ERROR", 1)) {
    return error_flag;
  }
  total_heat_release = NV_Ith_S(user_data->quadrature_,0);

  results[0] = current_time;
  results[1] = max_dp_dt;
  results[2] = time_max_dp_dt;
  results[3] = total_heat_release;
  results[4] = FindMaximum(time_history.size(),
                           &time_history[0],
                           1, // stride
                           &heat_release_rate_history[0],
                           1, // stride
                           false, // use quadratic interpolation
                           &time_at_history_max);
  results[5] = FindMaximum(time_history.size(),
                           &time_history[0],
                           1,
                           &heat_release_rate_history[0],
                           1,
                           true, // use quadratic interpolation
                           &time_at_history_max);
  ode_time_period = current_time - user_data->GetInitialTime();
  for(int j=1; j<user_data->num_quadratures_; ++j) {
    results[NUM_FIXED_RESULTS+j-1] =
      pow(NV_Ith_S(user_data->quadrature_,j)/ode_time_period,
          1.0/user_data->heat_release_rate_powers_[j]);
  }

  for(int j=0; j<num_burn_fractions; ++j) {
    results[j+NUM_FIXED_RESULTS+user_data->num_quadratures_-1] =
      burn_fraction_time[j];
  }

  (*simulation_time) = getHighResolutionTime()-sim_start;

  // write up the final stats
  if(thist_params.echo_stdout) {
    printf("# Final ODE time              [s]: %14.7e\n",current_time);
    printf("# Time at max dP/dt           [s]: %14.7e\n",time_max_dp_dt);
    printf("# Max dP/dt                [Pa/s]: %14.7e\n",max_dp_dt);
    printf("# Total Chemical Heat Release [J]: %14.7e\n",total_heat_release);
    printf("# Max heat release rate     [J/s]: %14.7e (linear)\n",
           results[4]);
    printf("# Max heat release rate     [J/s]: %14.7e (quadratic)\n",
           results[5]);
    printf("# Simulation wall clock time  [s]: %14.7e\n",*simulation_time);
    fflush(stdout);
    // if(thist != NULL) {
    //   fprintf(thist,"# Final ODE time             [s]: %14.7e\n",
    //           current_time);
    //   fprintf(thist,"# Time at max dP/dt          [s]: %14.7e\n",
    //           time_max_dp_dt);
    //   fprintf(thist,"# Max dP/dt               [Pa/s]: %14.7e\n",
    //           max_dp_dt);
    //   fprintf(thist,"# Simulation wall clock time [s]: %14.7e\n",
    //           *simulation_time);
    //   fflush(thist);
    // }
  }
  delete [] exceeded_burn_fraction;
  delete [] burn_fraction_time;
  return 0;
}


void BuildTimeHistoryLine(void *cvode_memory,
                          const double t,
                          N_Vector state,
                          UserData *user_data,
                          std::string &line)
{
  char format_line[MAX_LINE_LENGTH];
  double volume, volume_rate, temperature, pressure, density, total_moles, heat_release_rate;
  long int cv_num_steps, cv_num_rhs_evals, cv_num_lin_solve_setups;
  long int cv_num_err_test_fails, cv_dense_jac;
  int cv_last_order;
  double cv_last_step;
  std::vector<double> mole_fractions;

  line.clear();

  // state thermodynamics
  volume      = user_data->GetVolume(t, NV_DATA_S(state));
  volume_rate = user_data->GetVolumeRate(t, NV_DATA_S(state));
  temperature = user_data->GetTemperature(NV_DATA_S(state));
  pressure    = user_data->GetPressure(t,NV_DATA_S(state));
  density     = 0.0;
  total_moles = user_data->GetTotalMoles(NV_DATA_S(state));

  user_data->GetTrackedMoleFractions(NV_DATA_S(state),
                                    mole_fractions);

  heat_release_rate = ChemicalHeatReleaseRate(t,
                                              NV_DATA_S(state),
                                              user_data);


  sprintf(format_line,
          "%12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e",
          t,
          volume,
          volume_rate,
          temperature,
          pressure,
          total_moles,
          heat_release_rate);
  line += std::string(format_line);

  // cvode statistics
  cv_num_rhs_evals = 0;
  cv_num_lin_solve_setups = 0;
  cv_num_err_test_fails = 0;
  cv_last_order = 0;
  cv_last_step = 0.0;
  cv_dense_jac =0;
  CVodeGetNumSteps(cvode_memory,&cv_num_steps);
  if(cv_num_steps > 0) {
    CVodeGetNumRhsEvals(cvode_memory,&cv_num_rhs_evals);
    CVodeGetNumLinSolvSetups(cvode_memory,&cv_num_lin_solve_setups);
    CVodeGetNumErrTestFails(cvode_memory,&cv_num_err_test_fails);
    CVodeGetLastOrder(cvode_memory,&cv_last_order);
    CVodeGetLastStep(cvode_memory,&cv_last_step);
    // TODO: allow check for sparse and dense systems
#ifdef SUNDIALS4
    CVodeGetNumJacEvals(cvode_memory,&cv_dense_jac);
#else
    CVDlsGetNumJacEvals(cvode_memory,&cv_dense_jac);
#endif
  }
  sprintf(format_line,
	  "  %6d  %6d  %6d  %6d  %6d  %12.5e  %6d  %6d  %6d  %6d",
	  static_cast<int>(cv_num_steps),
          static_cast<int>(cv_num_rhs_evals),
          static_cast<int>(cv_num_lin_solve_setups),
          static_cast<int>(cv_num_err_test_fails),
          static_cast<int>(cv_dense_jac),
          cv_last_step,
	  static_cast<int>(cv_last_order),
	  user_data->num_preconditioner_setups,
	  user_data->num_new_jacobians,
	  user_data->num_preconditioner_solves);
  line += std::string(format_line);


  // mole fractions of tracked species
  for(size_t j=0; j<mole_fractions.size(); ++j) {
    sprintf(format_line,"  %12.5e",mole_fractions[j]);
    line += std::string(format_line);
  }

  line += std::string("\n");
}

void BuildTimeHistoryHeader(void *cvode_memory,
                            const double t,
                            N_Vector state,
                            UserData *user_data,
                            std::string &header)
{
  char format_line[MAX_LINE_LENGTH];
  std::vector<int> species_ids;
  user_data->GetTrackedSpeciesIds(species_ids);

  header.clear();
  header  = std::string("# Column  1: [s]     time\n");
  header += std::string("# Column  2: [m^3]   volume\n");
  header += std::string("# Column  3: [m^3/s] rate of volume change\n");
  header += std::string("# Column  4: [K]     temperature\n");
  header += std::string("# Column  5: [Pa]    pressure\n");
  header += std::string("# Column  6: [kmol]  total moles\n");
  header += std::string("# Column  7: [J/s]   heat release rate\n");
  header += std::string("# Column  8: [#]     time steps\n");
  header += std::string("# Column  9: [#]     rhs evaluations\n");
  header += std::string("# Column 10: [#]     linear solve setups\n");
  header += std::string("# Column 11: [#]     error test failures\n");
  header += std::string("# Column 12: [#]     Jacobian evaluations\n");
  header += std::string("# Column 13: [s]     last internal time step\n");
  header += std::string("# Column 14: [#]     last method order\n");
  header += std::string("# Column 15: [#]     number of preconditioner setups\n");
  header += std::string("# Column 16: [#]     number of new jacobians\n");
  header += std::string("# Column 17: [#]     number of preconditioner solves\n");

  // TODO add tracked species names
  for(size_t j=0; j<species_ids.size(); ++j) {
    sprintf(format_line,
            "# Column %2d: [-]     mole fraction based on %s (species id = %d)\n",
	    static_cast<int>(j+18),
            user_data->GetReactor()->GetNameOfStateId(species_ids[j]),
            species_ids[j]);
    header += std::string(format_line);
  }
}
void PrintStateDerivatives(const double t,
                           N_Vector state,
                           UserData *params)
{
  const int num_states = params->GetNumStates();
  N_Vector cvode_derivative;
  cvode_derivative = N_VNew_Serial(num_states);

  VariableVolumeRHS(t,state,cvode_derivative,params);
  printf("# state[id]             name: <state_value>  <state_derivative>\n");
  for(int j=0; j<num_states; ++j) {
    printf("# state[%2d] %16s: %.18g  %.18g\n",
           j,
           params->GetReactor()->GetNameOfStateId(j),
           NV_Ith_S(state,j),
           NV_Ith_S(cvode_derivative,j));
    fflush(stdout);
  }

  N_VDestroy_Serial(cvode_derivative);
}

int CheckFlag(void *flagvalue, const char *funcname, int opt)
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

  return 0;
}

static double FindMaximum(const size_t num_points,
                          const double x[],
                          const size_t x_stride,
                          const double f[],
                          const size_t f_stride,
                          const bool use_quadratic,
                          double *x_at_max)


{
  double interval_max, interval_x_at_max;
  double max_f = f[0];
  *x_at_max = x[0];
  if(use_quadratic && num_points >= 3) {
    for(size_t j=0; j<num_points-2; ++j) {
      interval_max = QuadraticMaximum(&x[j*x_stride],
                                      x_stride,
                                      &f[j*f_stride],
                                      f_stride,
                                      &interval_x_at_max);
      if(interval_max > max_f) {
        max_f = interval_max;
        *x_at_max = interval_x_at_max;
      }

    }
  }
  else {
    // max value at a grid point
    size_t x_id = x_stride;
    size_t f_id = f_stride;
    for(size_t j=1; j<num_points; ++j) {
      if(f[f_id] > max_f) {
        max_f = f[f_id];
        *x_at_max = x[x_id];
      }
      x_id += x_stride;
      f_id += f_stride;
    }
  }
  return max_f;
}

static double QuadraticMaximum(const double x[],
                               const size_t x_stride,
                               const double f[],
                               const size_t f_stride,
                               double *x_at_max)
{
  double x0 = x[0];
  double x1 = x[x_stride];
  double x2 = x[2*x_stride];
  double phi0 = f[0]         /((x0-x1)*(x0-x2));
  double phi1 = f[f_stride]  /((x1-x0)*(x1-x2));
  double phi2 = f[2*f_stride]/((x2-x0)*(x2-x1));

  double x_ext = 0.5*(phi0*(x1+x2) + phi1*(x0+x2) + phi2*(x0+x1))/
    (phi0+phi1+phi2);
  double f_ext;

  if(x0 < x_ext && x_ext < x2) {
    // extrema within domain may be maximum
    f_ext = phi0*(x_ext-x1)*(x_ext-x2) +
            phi1*(x_ext-x0)*(x_ext-x2) +
            phi2*(x_ext-x0)*(x_ext-x1);
    if(f_ext >= f[0] && f_ext >= f[f_stride] && f_ext >= f[2*f_stride]) {
      // should actually only need to compare to one value
      *x_at_max = x_ext;
      return f_ext;
    }
  } // otherwise maximum value is at the boundaries of the domain
  if(f[2*f_stride] > f[0]) {
    // last point
    *x_at_max = x[2*x_stride];
    return f[2*f_stride];
  }
  // first point
  *x_at_max = x[0];
  return f[0];
}
