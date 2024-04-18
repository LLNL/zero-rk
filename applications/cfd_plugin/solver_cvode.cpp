
#include <vector>

#include "solver_cvode.h"

#include "utility_funcs.h"

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#ifdef SUNDIALS2
#include <cvode/cvode_dense.h>
#include <cvode/cvode_spgmr.h>      // prototypes & constants for CVSPGMR
#elif SUNDIALS3
#include <cvode/cvode_direct.h>
#include <cvode/cvode_spils.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#elif SUNDIALS4
#include <sunmatrix/sunmatrix_dense.h>
#ifdef ZERORK_HAVE_SUNDIALS_LAPACK
#include <sunlinsol/sunlinsol_lapackdense.h>
#define ZERORK_USE_LAPACK_DENSE_MATRIX
#ifdef ZERORK_USE_LAPACK_DENSE_MATRIX
#include "interfaces/zrkmatrix/zrkmatrix_lapackdense.h"
#endif
#else
#include <sunlinsol/sunlinsol_dense.h>
#endif
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif

#ifdef ZERORK_GPU
  #ifdef ZERORK_HAVE_MAGMA
    #include <sunmatrix/sunmatrix_magmadense.h>
    #include <sunlinsol/sunlinsol_magmadense.h>
    #include <sunmemory/sunmemory_cuda.h>
    #include "interfaces/zrklinsol/zrklinsol_magmadense.h"
  #endif
#endif

CvodeSolver::CvodeSolver(ReactorBase& reactor)
  :
      SolverBase(reactor),
      reactor_ref_(reactor),
      cb_fn_(nullptr),
      cb_fn_data_(nullptr)
{}

int CvodeSolver::Integrate(const double end_time) {
  N_Vector& state = reactor_ref_.GetStateNVectorRef();
  N_Vector derivative = N_VClone(state);
  int reactor_id = reactor_ref_.GetID();

#if defined SUNDIALS3 || defined SUNDIALS4
  SUNMatrix A;
  SUNLinearSolver LS;
#endif
#ifdef SUNDIALS4
  SUNNonlinearSolver NLS;
  void* cvode_mem = CVodeCreate(CV_BDF);
#else
  void* cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#endif
  int flag = CVodeInit(cvode_mem, ReactorGetTimeDerivative, 0.0, state);
  check_cvode_flag(&flag, "CVodeInit", 1);

  flag = CVodeSetUserData(cvode_mem, &(reactor_ref_));
  check_cvode_flag(&flag, "CVodeSetUserData", 1);

  N_Vector abs_tol_vector = N_VClone(state);
  N_Vector abs_tol_corrections = N_VClone(state);
  reactor_ref_.GetAbsoluteToleranceCorrection(abs_tol_corrections);
  N_VConst(double_options_["abs_tol"],abs_tol_vector);
  N_VProd(abs_tol_vector, abs_tol_corrections, abs_tol_vector);
  flag = CVodeSVtolerances(cvode_mem, double_options_["rel_tol"], abs_tol_vector);
  check_cvode_flag(&flag, "CVodeSVtolerances", 1);
  N_VDestroy(abs_tol_corrections);

  flag = CVodeSetMaxNumSteps(cvode_mem, int_options_["max_steps"]);
  check_cvode_flag(&flag, "CVodeSetMaxNumSteps", 1);

  flag = CVodeSetMaxStep(cvode_mem, double_options_["max_dt"]);
  check_cvode_flag(&flag, "CVodeSetMaxStep", 1);

  flag = CVodeSetNonlinConvCoef(cvode_mem, double_options_["nonlinear_convergence_coeff"]);
  check_cvode_flag(&flag, "CVodeSetNonlinConvCoef", 1);

  if(reactor_ref_.GetNumRootFunctions() > 0) {
      flag = CVodeRootInit(cvode_mem, reactor_ref_.GetNumRootFunctions(), ReactorRootFunction);
      check_cvode_flag(&flag, "CVodeRootInit", 1);
  }

  int num_batches = reactor_ref_.GetNumBatchReactors();
  reactor_ref_.GetReactorWeightsRef().assign(num_batches,1.0);

  bool dense_direct = (int_options_["dense"]==1 && int_options_["iterative"]==0);
  if(dense_direct) {
#ifdef ZERORK_GPU
  #ifdef ZERORK_HAVE_MAGMA
    //Check for GPU
    if(num_batches > 1) {
        int NEQ = reactor_ref_.GetNumStateVariables();
        SUNMemoryHelper memhelper = SUNMemoryHelper_Cuda();
        //NULL is passed to use default cuda stream
        A = SUNMatrix_MagmaDenseBlock(num_batches, NEQ, NEQ, SUNMEMTYPE_DEVICE, memhelper, NULL);
        LS = ZRKLinSol_MagmaDense(state, A);
        ZRKLinSol_MagmaDense_SetAsync(LS, 0);
        flag = CVodeSetLinearSolver(cvode_mem, LS, A);
        check_cvode_flag(&flag, "CVodeSetLinearSolver", 1);
        flag = CVodeSetJacFn(cvode_mem, ReactorGetJacobianDense);
        check_cvode_flag(&flag, "CVodeSetJacFn", 1);
    } else
  #endif
#endif
    {
#ifdef SUNDIALS2
      flag = CVDense(cvode_mem, reactor_ref_.GetNumStateVariables());
      check_cvode_flag(&flag, "CVLapackDense", 1);
      if(int_options_["analytic"]) {
        flag = CVDlsSetDenseJacFn(cvode_mem, ReactorGetJacobianDense);
        check_cvode_flag(&flag, "CVDlsDenseJacFn", 1);
      } else {
        flag = CVDlsSetDenseJacFn(cvode_mem, NULL);
        check_cvode_flag(&flag, "CVDlsDenseJacFn", 1);
      }
#elif defined SUNDIALS3
      int NEQ = reactor_ref_.GetNumStateVariables();
      A = SUNDenseMatrix(NEQ, NEQ);
      LS = SUNDenseLinearSolver(state, A);
      flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
      check_cvode_flag(&flag, "CVDlsSetLinearSolver", 1);
      if(int_options_["analytic"]) {
        flag = CVDlsSetJacFn(cvode_mem, ReactorGetJacobianDense);
        check_cvode_flag(&flag, "CVDlsSetJacFn", 1);
      } else {
        flag = CVDlsSetJacFn(cvode_mem, NULL);
        check_cvode_flag(&flag, "CVDlsSetJacFn", 1);
      }
#elif defined SUNDIALS4
      int NEQ = reactor_ref_.GetNumStateVariables();
#ifdef ZERORK_HAVE_SUNDIALS_LAPACK
#ifdef ZERORK_USE_LAPACK_DENSE_MATRIX
      A = ZRKDenseLapackMatrix(NEQ, NEQ);
#else
      A = SUNDenseMatrix(NEQ, NEQ);
#endif
      LS = SUNLinSol_LapackDense(state, A);
#else
      A = SUNDenseMatrix(NEQ, NEQ);
      LS = SUNLinSol_Dense(state, A);
#endif
      flag = CVodeSetLinearSolver(cvode_mem, LS, A);
      check_cvode_flag(&flag, "CVodeSetLinearSolver", 1);
      if(int_options_["analytic"]) {
        flag = CVodeSetJacFn(cvode_mem, ReactorGetJacobianDense);
        check_cvode_flag(&flag, "CVodeSetJacFn", 1);
      } else {
        flag = CVodeSetJacFn(cvode_mem, NULL);
        check_cvode_flag(&flag, "CVodeSetJacFn", 1);
      }
#else
#error "Unknown SUNDIALS version"
#endif
    }
  } else {
#ifdef SUNDIALS2
    flag = CVSpgmr(cvode_mem, PREC_LEFT, 0);
    check_cvode_flag(&flag, "CVSpgmr", 1);

    flag = CVSpilsSetGSType(cvode_mem, MODIFIED_GS);
    check_cvode_flag(&flag, "CVSpilsSetGSType", 1);

    flag = CVSpilsSetPreconditioner(cvode_mem,
                                    ReactorJacobianSetupAndFactor,
                                    ReactorJacobianSolveCVODE);
    check_cvode_flag(&flag, "CVSpilsSetPreconditioner", 1);
#elif defined SUNDIALS3
    LS = SUNSPGMR(state, PREC_LEFT, 0);

    flag = SUNSPGMRSetGSType(LS, MODIFIED_GS);
    check_cvode_flag(&flag, "SUNSPGRMSetGSType", 1);

    flag = CVSpilsSetLinearSolver(cvode_mem, LS);
    check_cvode_flag(&flag, "CVSpilsSetLinearSolver", 1);

    flag = CVSpilsSetPreconditioner(cvode_mem,
                                    ReactorJacobianSetupAndFactor,
                                    ReactorJacobianSolveCVODE);
    check_cvode_flag(&flag, "CVSpilsSetPreconditioner", 1);
#elif defined SUNDIALS4
    NLS = SUNNonlinSol_Newton(state);
    flag = CVodeSetNonlinearSolver(cvode_mem, NLS);
    check_cvode_flag(&flag, "CVodeSetNonlinearSolver", 1);

    LS = SUNLinSol_SPGMR(state, PREC_LEFT, 0);
    flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    check_cvode_flag(&flag, "CVodeSetLinearSolver", 1);

    flag = CVodeSetPreconditioner(cvode_mem,
                                  ReactorJacobianSetupAndFactor,
                                  ReactorJacobianSolveCVODE);
    check_cvode_flag(&flag, "CVodeSetPreconditioner", 1);
#else
#error "Unknown SUNDIALS version"
#endif

#if defined SUNDIALS2 || defined SUNDIALS3
    flag = CVSpilsSetEpsLin(cvode_mem, double_options_["eps_lin"]);    // Default [0.05]
#else
    flag = CVodeSetEpsLin(cvode_mem, double_options_["eps_lin"]);    // Default [0.05]
#endif
    check_cvode_flag(&flag, "CVSpilsSetEpsLin", 1);
  }

  long int nsteps = 0;
  long int nsteps_curr = 0;
  long int num_linear_solve_setups = 0;
  const int max_tries = int_options_["cvode_num_retries"];
  const double abs_tol_adj = double_options_["cvode_retry_absolute_tolerance_adjustment"];
  const double rel_tol_adj = double_options_["cvode_retry_relative_tolerance_adjustment"];
  double rel_tol = double_options_["rel_tol"];
  int num_tries = 0;
  double tcurr = 0.0;
  double tprev = 0.0;
  int cb_flag = 0;
  int cv_mode = CV_ONE_STEP;
  while(tcurr < end_time) {
      flag = CVode(cvode_mem, end_time, state, &tcurr, cv_mode);
      CVodeGetNumSteps(cvode_mem,&nsteps_curr);
      nsteps += 1;
      if(nsteps_curr >= int_options_["max_steps"]) {
        flag = CV_TOO_MUCH_WORK;
      }
      if(check_cvode_flag(&flag, "CVode", 1)) {
        num_tries += 1;
        N_VScale(abs_tol_adj, abs_tol_vector, abs_tol_vector);
        rel_tol *= rel_tol_adj;
        CVodeSVtolerances(cvode_mem, rel_tol, abs_tol_vector);

        if(num_tries == max_tries) {
          break;
        }
        CVodeReInit(cvode_mem, tcurr, state);
      }
      long int last_nlss = num_linear_solve_setups;
      CVodeGetNumLinSolvSetups(cvode_mem, &num_linear_solve_setups);
      if( (flag == CV_SUCCESS || flag == CV_ROOT_RETURN) && cb_fn_ != nullptr ) {
        if(N_VGetVectorID(state) == SUNDIALS_NVEC_SERIAL) {
          flag = CVodeGetDky(cvode_mem, tcurr, 1, derivative);
          double dt = tcurr - tprev;
          cb_flag = cb_fn_(reactor_id, nsteps, tcurr, dt, NV_DATA_S(state), NV_DATA_S(derivative), cb_fn_data_);
          if(cb_flag != 0) {
            break;
          }
        }
      }
      if(flag == CV_ROOT_RETURN) {
        flag = CV_SUCCESS;
        reactor_ref_.SetRootTime(tcurr);
        if(int_options_["stop_after_ignition"]) {
          break;
        }
      }
      tprev = tcurr;
  }
  if(cv_mode == CV_ONE_STEP && flag == CV_SUCCESS) {
    flag = CVodeGetDky(cvode_mem, std::min(tcurr,end_time), 0, state);
  }

  if(flag < 0) {
    printf("WARNING: Failed to complete integration.\n");
    if(nsteps <= 0) {
      nsteps = -1;
    } else {
      nsteps = -nsteps;
    }
  }

  //Cleanup memory
#if defined SUNDIALS3 || defined SUNDIALS4
  SUNLinSolFree(LS);
  if(dense_direct) {
    SUNMatDestroy(A);
  }
#endif
#if defined SUNDIALS4
  if(!dense_direct) {
    SUNNonlinSolFree(NLS);
  }
#endif
  N_VDestroy(derivative);
  N_VDestroy(abs_tol_vector);

  CVodeFree(&cvode_mem);

  return nsteps;
}

int CvodeSolver::Iterative() {
  if(int_options_["dense"] == 1) {
    return int_options_["iterative"];
  } else {
    return 1;
  }
}

void CvodeSolver::AdjustWeights(void* cvode_mem) {
  int num_batches = reactor_ref_.GetNumBatchReactors();
  std::vector<double>& weights = reactor_ref_.GetReactorWeightsRef();
  if(num_batches == 1) {
    weights[0] = 1.0;
    return;
  }

  //Get error weights and local errors.
  N_Vector& state = reactor_ref_.GetStateNVectorRef();
  N_Vector eweight = N_VClone(state);
  N_Vector ele = N_VClone(state);
  N_Vector error = N_VClone(state);
  N_Vector error_batch = N_VClone(state);
  int flag = CVodeGetErrWeights(cvode_mem, eweight);
  flag = CVodeGetEstLocalErrors(cvode_mem, ele); 
  N_VProd(eweight, ele, error);
  //N_VAbs(error, error);

  std::vector<double> max_errors(num_batches,0.0);
  N_Vector batch_mask = N_VClone(state);
  for(int k = 0; k < num_batches; ++k) {
    reactor_ref_.SetBatchMaskNVector(k,batch_mask);
    N_VProd(batch_mask, error, error_batch);
    max_errors[k] = N_VMaxNorm(error_batch);
  }
  N_VDestroy(batch_mask);
  std::vector<double> sorted_max_errors = max_errors;
  std::sort(sorted_max_errors.begin(),sorted_max_errors.end());
  //double min_max_err = sorted_max_errs[0];
  double max_max_err = sorted_max_errors[num_batches-1];
  double med_max_err = sorted_max_errors[num_batches/2-1];
  double avg_max_err = std::accumulate(max_errors.begin(),max_errors.end(),0.0)/num_batches;

  //TODO: More robust factors for error multipliers
  for(int k = 0; k < num_batches; ++k) {
      if(max_errors[k] == max_max_err) {
          weights[k] *= 1.00;
      } else if(max_errors[k] > avg_max_err) {
          weights[k] *= 0.9999;
      } else if(max_errors[k] > med_max_err) {
          weights[k] *= 0.999;
      } else {
          weights[k] *= 0.99;
      }
  }
  N_VDestroy(eweight);
  N_VDestroy(ele);
  N_VDestroy(error);
  N_VDestroy(error_batch);
}

void CvodeSolver::SetCallbackFunction(zerork_callback_fn fn, void* cb_fn_data) {
  cb_fn_ = fn;
  cb_fn_data_ = cb_fn_data;
}

