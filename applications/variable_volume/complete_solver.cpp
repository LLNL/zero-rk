#include <stdlib.h>
#include <stdio.h>

#include <string>
#include <vector>

#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros

#ifdef SUNDIALS2
#include <cvode/cvode_dense.h>      // prototypes for CVDense
#include <cvode/cvode_spgmr.h>      // prototypes & constants for CVSPGMR
#elif SUNDIALS3
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <cvode/cvode_spils.h>      // prototypes & constants for CVSPGMR
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#elif SUNDIALS4
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif

#include "user_functions.h"
#include "complete_solver.h"

const int MAX_LINE_LENGTH = 8192;

using zerork::getHighResolutionTime;

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
  if(user_data->GetParser()->negativeConcentrationMitigation()) {
    error_flag = CVodeInit(cvode_memory,
                           VariableVolumeRHS_Limit1,
                           user_data->GetInitialTime(),
                           cvode_state);
  } else {
    error_flag = CVodeInit(cvode_memory,
                           VariableVolumeRHS,
                           user_data->GetInitialTime(),
                           cvode_state);
  }
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

  error_flag = CVodeSetStabLimDet(cvode_memory,
                                  TRUE);
  if (CheckFlag(&error_flag, "CVodeSetStabLimDet", 1)) exit(-1);

  // Attach the linear solver
  if(user_data->GetParser()->Jacobian() == std::string("DENSE_ANALYTIC") ||
     user_data->GetParser()->Jacobian() == std::string("DENSE_NUMERIC")) {
    SetupDenseSolver(cvode_memory,
                     aux_cvode,
                     user_data,
                     cvode_state);
  } else {
    SetupSparseSolver(cvode_memory,
                      aux_cvode,
                      user_data,
                      cvode_state);
  }
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
                 5); // default is five
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
#ifdef SUNDIALS2
  // Dense solver
  error_flag = CVDense(cvode_memory,
                       user_data->GetNumStates());
  if (CheckFlag(&error_flag, "CVDense", 1)) exit(-1);

  // Attach the dense jacobian function
  if(user_data->GetParser()->Jacobian() == std::string("DENSE_ANALYTIC")) {
    error_flag = CVDlsSetDenseJacFn(cvode_memory,
                                    VariableVolumeDenseJacobian);
    if (CheckFlag(&error_flag, "CVDlsSetDenseJacFn", 1)) exit(-1);
  }
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

  // ------------------------------------------------------------------------
  // set the initial time and state, and the final time
  //user_data.GetInitialState(NV_DATA_S(cvode_state));
  current_time = user_data->GetInitialTime();
  final_time = user_data->GetParser()-> finalTime();

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

  max_dp_dt      = 0.0;
  time_max_dp_dt = current_time;
  current_pressure =
    user_data->GetPressure(current_time,NV_DATA_S(cvode_state));

  // set the last print more than one min_time_step in the past
  print_time = current_time - 1.1*thist_params.min_time_step;


  while(current_time < final_time) {
    previous_time = current_time;
    previous_pressure = current_pressure;

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
          fflush(stdout);
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
  }
  results[0] = current_time;
  results[1] = max_dp_dt;
  results[2] = time_max_dp_dt;
  (*simulation_time) = getHighResolutionTime()-sim_start;

  // write up the final stats
  if(thist_params.echo_stdout) {
    printf("# Final ODE time             [s]: %14.7e\n",current_time);
    printf("# Time at max dP/dt          [s]: %14.7e\n",time_max_dp_dt);
    printf("# Max dP/dt               [Pa/s]: %14.7e\n",max_dp_dt);
    printf("# Simulation wall clock time [s]: %14.7e\n",*simulation_time);
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

  return 0;
}


void BuildTimeHistoryLine(void *cvode_memory,
                          const double t,
                          N_Vector state,
                          UserData *user_data,
                          std::string &line)
{
  char format_line[MAX_LINE_LENGTH];
  double volume, volume_rate, temperature, pressure, total_moles;
  double heat_release_rate;
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
  total_moles = user_data->GetTotalMoles(NV_DATA_S(state));
  heat_release_rate = user_data->GetChemicalHeatReleaseRate(t,NV_DATA_S(state));

  user_data->GetTrackedMoleFractions(NV_DATA_S(state),
                                    mole_fractions);

  const char * format_str;
  if(user_data->GetParser()->timeHistoryOutputHighPrecision()) {
    format_str = "%24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e";
  } else {
    format_str = "%12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e";
  }

  sprintf(format_line,
          format_str,
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
    //CVDlsGetNumJacEvals(cvode_memory,&cv_dense_jac);
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

  // track the minimum mole fraction
  const int num_species = user_data->GetReactor()->GetNumStates()-2;
  mole_fractions.clear();
  mole_fractions.assign(num_species, 0.0);
  user_data->GetMoleFractions(NV_DATA_S(state),
                              &mole_fractions[0]);
  double min_mole_fraction = mole_fractions[0];
  int min_species_id = 0;
  for(int j=1; j<num_species; ++j) {
    if(mole_fractions[j] < min_mole_fraction) {
      min_mole_fraction = mole_fractions[j];
      min_species_id = j;
    }
  }
  sprintf(format_line,"  %12.5e  (%s id=%d)",
          min_mole_fraction,
          user_data->GetReactor()->GetNameOfStateId(min_species_id),
          min_species_id);
  line += std::string(format_line);


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
  header += std::string("# Column  7: [J/s]   chemical heat release rate\n");
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
  sprintf(format_line,
          "# Column %2d: [-]     minimum mole fraction\n",
	    static_cast<int>(species_ids.size()+18));
    header += std::string(format_line);

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
