#ifndef SOLVER_BASE_H_
#define SOLVER_BASE_H_

#include "optionable.h"
#include "reactor_base.h"

class SolverBase : public Optionable
{
 public:
  SolverBase(ReactorBase& reactor) {};
  virtual ~SolverBase() {};

  virtual int Integrate(const double time) = 0;
  virtual int Iterative() = 0;
};

#endif
