#ifndef SOLVER_SODEX_H_
#define SOLVER_SODEX_H_

#include <string>

#include "reactor_base.h"
#include "solver_base.h"
#include "sundials/sundials_nvector.h"

class SodexSolver : public SolverBase {
 public:
  SodexSolver(ReactorBase& reactor);
  ~SodexSolver(){};

  int Integrate(const double end_time);
  int Iterative() { return 0; };

  void SetCallbackFunction(zerork_callback_fn fn, void* cb_fn_data);
  int MonitorFn(int nsteps, double x, double h, N_Vector y, N_Vector ydot);

 private:
  ReactorBase& reactor_ref_;
  void AdjustWeights();

  zerork_callback_fn cb_fn_;
  void* cb_fn_data_;
};

#endif
