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

 private:
  ReactorBase& reactor_ref_;
  void AdjustWeights();
  static int RootMonitor(int nsteps, double x, N_Vector y, N_Vector ydot, void *user_data);
};

#endif
