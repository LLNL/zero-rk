#ifndef ZERORK_REACTOR_MANAGER_BASE_H
#define ZERORK_REACTOR_MANAGER_BASE_H

#include <string>

#include "optionable.h"

#include "zerork_cfd_plugin.h" // for zerork_field_t, zerork_status_t

class ZeroRKReactorManagerBase : public Optionable
{
 public:
  ZeroRKReactorManagerBase() {};
  virtual ~ZeroRKReactorManagerBase() {};

  virtual zerork_status_t ReadOptionsFile(const std::string& options_filename) = 0;
  virtual zerork_status_t LoadMechanism() = 0;

  virtual zerork_status_t SetInputVariables(int n_cycle,
                                double time,
                                double dt,
                                int n_reactors,
                                double* T,
                                double* P,
                                double* mass_fractions) = 0;

  virtual zerork_status_t SetAuxFieldPointer(zerork_field_t ft, double* field_pointer) = 0;
  virtual zerork_status_t SetCallbackFunction(zerork_callback_fn fn, void* user_data) = 0;
  virtual zerork_status_t SetReactorIDs(int* field_pointer) = 0;

  virtual zerork_status_t FinishInit() = 0;
  virtual zerork_status_t LoadBalance() = 0;
  virtual zerork_status_t SolveReactors() = 0;
  virtual zerork_status_t RedistributeResults() = 0;
  virtual zerork_status_t PostSolve() = 0;
};

#endif
