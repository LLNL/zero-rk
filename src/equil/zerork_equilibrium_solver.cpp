

#include "zerork_equilibrium_solver.h"

namespace zerork {

equilibrium_solver::equilibrium_solver(const zerork::mechanism& mech)
#ifdef ZERORK_HAVE_EQUILIBRIUM_SOLVER
 : zcm_(mech)
#endif
{};


void equilibrium_solver::equilibrate(double p, double& T, double* y) {
#ifdef ZERORK_HAVE_EQUILIBRIUM_SOLVER
  zcm_.equilibrate(p, T, y);
#else
  throw std::runtime_error("No equilbrium solver available.");
#endif
}

equilibrium_solver::~equilibrium_solver() {};

} //namespace zerork

