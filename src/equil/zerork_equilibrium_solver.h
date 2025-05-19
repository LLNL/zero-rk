
#ifndef ZERORK_EQUILIBRIUM_SOLVER_H
#define ZERORK_EQUILIBRIUM_SOLVER_H

#include "mechanism.h"

#undef ZERORK_HAVE_EQUILIBRIUM_SOLVER
#ifdef ZERORK_HAVE_CEQ
#include "interfaces/ceq/zerork_ceq_manager.h"
#define ZERORK_HAVE_EQUILIBRIUM_SOLVER
#endif
#ifdef ZERORK_HAVE_CANTERA
#include "interfaces/cantera/zerork_cantera_manager.h"
#define ZERORK_HAVE_EQUILIBRIUM_SOLVER
#endif

namespace zerork {

class equilibrium_solver {

 public:

  equilibrium_solver(const zerork::mechanism& mech);
  ~equilibrium_solver();

  void equilibrate(double p, double& T, double* y);

 private:
#ifdef ZERORK_HAVE_CEQ
  zerork::ceq_manager zcm_;
#else
#ifdef ZERORK_HAVE_CANTERA
  zerork::cantera_manager zcm_;
#endif
#endif
};

} //namespace zerork

#endif
