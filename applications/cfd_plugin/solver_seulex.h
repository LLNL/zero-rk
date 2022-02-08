#ifndef SOLVER_SEULEX_H_
#define SOLVER_SEULEX_H_

#include <string>

#include "solver_base.h"
#include "reactor_base.h"

#include "sundials/sundials_nvector.h"

class SeulexSolver : public SolverBase
{
 public:
  SeulexSolver(ReactorBase& reactor);
  ~SeulexSolver() {};

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
