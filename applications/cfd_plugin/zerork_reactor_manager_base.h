#ifndef ZERORK_REACTOR_MANAGER_BASE_H
#define ZERORK_REACTOR_MANAGER_BASE_H

#include <string>

#include "optionable.h"

#include "zerork_cfd_plugin.h" // for zerork_field_type

class ZeroRKReactorManagerBase : public Optionable
{
 public:
  ZeroRKReactorManagerBase(const char *input_filename,
                           const char* mech_filename,
                           const char* therm_filename) {};
  virtual ~ZeroRKReactorManagerBase() {};

  virtual void SetInputVariables(int n_cycle,
                                double time,
                                double dt,
                                int n_reactors, 
                                double* T,
                                double* P,
                                double* mass_fractions) = 0;

  virtual void SetAuxFieldPointer(zerork_field_type ft, double* field_pointer) = 0;

  virtual void LoadBalance() = 0;
  virtual void SolveReactors() = 0;
  virtual void RedistributeResults() = 0;
  virtual void PostSolve() = 0; 
};

#endif
