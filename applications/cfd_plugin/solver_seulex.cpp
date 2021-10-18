
#include <vector>

#include "solver_seulex.h"
#include "utility_funcs.h"
#include "interfaces/seulex_cpp/seulex_cpp.h"

SeulexSolver::SeulexSolver(ReactorBase& reactor)
  :
      SolverBase(reactor),
      reactor_ref_(reactor)
{}

int SeulexSolver::Integrate(const double end_time) {
  N_Vector& state = reactor_ref_.GetStateNVectorRef();
  int num_variables = reactor_ref_.GetNumStateVariables();
  int num_batches = reactor_ref_.GetNumBatchReactors();
  seulex_cpp::seulex s(num_variables*num_batches, state);

  N_Vector abs_tol_vector = N_VClone(state);
  N_Vector rel_tol_vector = N_VClone(state);
  N_VConst(double_options_["abs_tol"],abs_tol_vector);
  N_VConst(double_options_["rel_tol"],rel_tol_vector);

  N_Vector abs_tol_corrections = N_VClone(state);
  reactor_ref_.GetAbsoluteToleranceCorrection(abs_tol_corrections);
  N_VProd(abs_tol_vector, abs_tol_corrections, abs_tol_vector);

  s.set_tolerances(rel_tol_vector,abs_tol_vector);
  N_VDestroy(abs_tol_vector);
  N_VDestroy(rel_tol_vector);
  N_VDestroy(abs_tol_corrections);

  s.set_hinit(end_time);
  s.set_nmax(int_options_["max_steps"]);
  s.set_hmax(std::min(end_time,double_options_["max_dt"]));
  s.set_user_data((void*) &(reactor_ref_));
  s.set_work_params(1.0, 3.0, 8.0, 1.5); //fcn, jac, decomp, solve
  s.set_nsequ(2);

  s.set_deriv_fcn(ReactorGetTimeDerivative);
  s.set_jac_fcn(ReactorJacobianSetup);
  s.set_jac_decomp_fcn(ReactorJacobianFactor);
  s.set_jac_solve_fcn(ReactorJacobianSolve);

  s.set_output_fcn(this->RootMonitor);

  reactor_ref_.GetReactorWeightsRef().assign(num_batches,1.0);

  int flag = 0;
  int maxRetry = 3;
  int nRetry   = 0;
  double tcurr = 0.0;
  while(tcurr < end_time) {
      flag = s.solve(&tcurr,end_time,state);
      if (flag != 0) {
        //printf("WARNING: Failed integration current temperature: %g\n",NV_Ith_S(systemState,systemParam.nSpc));
        if(nRetry < maxRetry) {
          ++nRetry;
          continue;
        } else {
          break;
        }
      }
      //if( last_nlss != num_linear_solve_setups ) {
      //if( nsteps % 100 == 0) {
        //AdjustWeights();
      //}
      //}
      if (flag == 0) break;
  }
  int nfcn, njac, nsteps, naccpt, nrejct, ndec, nsol;
  s.get_integrator_stats(&nfcn,&njac,&nsteps,&naccpt,&nrejct,&ndec,&nsol);

  if(flag != 0) {
    printf("WARNING: Failed to complete integration.\n");
    if(nsteps <= 0) {
      nsteps = -1;
    } else {
      nsteps = -nsteps;
    }
  }

  return nsteps;
}


void SeulexSolver::AdjustWeights() {
  int num_batches = reactor_ref_.GetNumBatchReactors();
  std::vector<double>& weights = reactor_ref_.GetReactorWeightsRef();
  if(num_batches == 1) {
    weights[0] = 1.0;
    return;
  }
  //Pass
  return;
}


int SeulexSolver::RootMonitor(int nsteps, double x, N_Vector y, N_Vector ydot, void *user_data)
{
  int num_root_fns = ReactorGetNumRootFunctions(user_data);
  static std::vector<double> last_root_fn_values(num_root_fns,0.0);
  std::vector<double> current_root_fn_values(num_root_fns);
  int flag = ReactorRootFunction(x, y, &current_root_fn_values[0], user_data);
  if(nsteps > 0) {
    for (int i = 0; i < num_root_fns; ++i) {
      if(current_root_fn_values[i] * last_root_fn_values[i] <= 0) {
          //printf("Found root[%d]: %g  (nsteps=%d)\n",i,x,nsteps);
      }
    }
  }
  for (int i = 0; i < num_root_fns; ++i) {
    last_root_fn_values[i] = current_root_fn_values[i];
  }
  return 0;
}

