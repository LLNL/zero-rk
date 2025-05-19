
#include "solver_sodex.h"

#include <vector>

#include "interfaces/sodex_cpp/sodex_cpp.h"
#include "nvector/nvector_serial.h"
#include "utility_funcs.h"

static int SodexSolverMonitorFn(int nsteps, double x, double h, N_Vector y,
                                N_Vector ydot, void* cb_fn_data) {
  SodexSolver* solver = static_cast<SodexSolver*>(cb_fn_data);
  return solver->MonitorFn(nsteps, x, h, y, ydot);
}

SodexSolver::SodexSolver(ReactorBase& reactor)
    : SolverBase(reactor),
      reactor_ref_(reactor),
      cb_fn_(nullptr),
      cb_fn_data_(nullptr) {}

int SodexSolver::Integrate(const double end_time) {
  N_Vector& state = reactor_ref_.GetStateNVectorRef();
  int num_variables = reactor_ref_.GetNumStateVariables();
  int num_batches = reactor_ref_.GetNumBatchReactors();
  sodex_cpp::sodex s(num_variables * num_batches, state);

  N_Vector abs_tol_vector = N_VClone(state);
  N_Vector rel_tol_vector = N_VClone(state);
  N_VConst(double_options_["abs_tol"], abs_tol_vector);
  N_VConst(double_options_["rel_tol"], rel_tol_vector);

  N_Vector abs_tol_corrections = N_VClone(state);
  reactor_ref_.GetAbsoluteToleranceCorrection(abs_tol_corrections);
  N_VProd(abs_tol_vector, abs_tol_corrections, abs_tol_vector);

  s.set_tolerances(rel_tol_vector, abs_tol_vector);
  N_VDestroy(abs_tol_vector);
  N_VDestroy(rel_tol_vector);
  N_VDestroy(abs_tol_corrections);

  s.set_hinit(end_time);
  s.set_nmax(int_options_["max_steps"]);
  s.set_hmax(std::min(end_time, double_options_["max_dt"]));
  s.set_user_data((void*)&(reactor_ref_));
  s.set_work_params(1.0, 3.0, 8.0, 1.5);  // fcn, jac, decomp, solve
  s.set_nsequ(1);

  s.set_deriv_fcn(ReactorGetTimeDerivative);
  s.set_jac_fcn(ReactorJacobianSetup);
  s.set_jac_decomp_fcn(ReactorJacobianFactor);
  s.set_jac_solve_fcn(ReactorJacobianSolve);

  s.set_output_fcn(SodexSolverMonitorFn, this);

  reactor_ref_.GetReactorWeightsRef().assign(num_batches, 1.0);

  int flag = 0;
  int maxRetry = 3;
  int nRetry = 0;
  double tcurr = 0.0;
  while (tcurr < end_time) {
    flag = s.solve(&tcurr, end_time, state);
    if (flag != 0) {
      if (nRetry < maxRetry) {
        ++nRetry;
        continue;
      } else {
        break;
      }
    }
    if (flag == 0) break;
  }
  int nfcn, njac, nsteps, naccpt, nrejct, ndec, nsol;
  s.get_integrator_stats(&nfcn, &njac, &nsteps, &naccpt, &nrejct, &ndec, &nsol);

  if (flag != 0) {
    printf("WARNING: Failed to complete integration.\n");
    if (nsteps <= 0) {
      nsteps = -1;
    } else {
      nsteps = -nsteps;
    }
  }

  return nsteps;
}

void SodexSolver::AdjustWeights() {
  int num_batches = reactor_ref_.GetNumBatchReactors();
  std::vector<double>& weights = reactor_ref_.GetReactorWeightsRef();
  if (num_batches == 1) {
    weights[0] = 1.0;
    return;
  }
  // Pass
  return;
}

void SodexSolver::SetCallbackFunction(zerork_callback_fn fn, void* cb_fn_data) {
  cb_fn_ = fn;
  cb_fn_data_ = cb_fn_data;
}

int SodexSolver::MonitorFn(int nsteps, double x, double h, N_Vector y,
                           N_Vector ydot) {
  int cb_flag = 0;
  int reactor_id = reactor_ref_.GetID();
  int num_root_fns = ReactorGetNumRootFunctions(&reactor_ref_);
  if (num_root_fns > 0) {
    static std::vector<double> last_root_fn_values(num_root_fns, 0.0);
    std::vector<double> current_root_fn_values(num_root_fns);
    int flag =
        ReactorRootFunction(x, y, &current_root_fn_values[0], &reactor_ref_);
    if (nsteps > 0) {
      for (int i = 0; i < num_root_fns; ++i) {
        if (current_root_fn_values[i] * last_root_fn_values[i] <= 0) {
          // printf("Found root[%d]: %g  (nsteps=%d)\n",i,x,nsteps);
        }
      }
    }
    for (int i = 0; i < num_root_fns; ++i) {
      last_root_fn_values[i] = current_root_fn_values[i];
    }
  }
  if (cb_fn_ != nullptr) {
    if (N_VGetVectorID(y) == SUNDIALS_NVEC_SERIAL) {
      cb_flag = cb_fn_(reactor_id, nsteps, x, h, NV_DATA_S(y), NV_DATA_S(ydot),
                       cb_fn_data_);
    }
  }
  return cb_flag;
}
