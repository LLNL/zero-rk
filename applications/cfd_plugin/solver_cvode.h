#ifndef SOLVER_CVODE_H_
#define SOLVER_CVODE_H_

#include <string>

#include "solver_base.h"
#include "reactor_base.h"

#include "sundials/sundials_nvector.h"

class CvodeSolver : public SolverBase
{
 public:
  CvodeSolver(ReactorBase& reactor);
  ~CvodeSolver() {};

  int Integrate(const double end_time);
  int Iterative();

 private:
  ReactorBase& reactor_ref_;
  void AdjustWeights(void* cvode_mem);
};

#endif
