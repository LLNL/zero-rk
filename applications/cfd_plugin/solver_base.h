#ifndef SOLVER_BASE_H_
#define SOLVER_BASE_H_

#include "optionable.h"
#include "reactor_base.h"

#include "zerork_cfd_plugin.h" //for zerork_callback_fn

class SolverBase : public Optionable
{
 public:
  SolverBase(ReactorBase& reactor) {};
  virtual ~SolverBase() {};

  virtual int Integrate(const double time) = 0;
  virtual int Iterative() = 0;

  virtual void SetCallbackFunction(zerork_callback_fn fn, void* user_data) = 0;
};

#endif
