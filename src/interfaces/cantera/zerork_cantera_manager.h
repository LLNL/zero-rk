
#ifndef ZERORK_CANTERA_MANAGER_H
#define ZERORK_CANTERA_MANAGER_H

#include "mechanism.h"
#include "cantera/thermo/IdealGasPhase.h"

namespace zerork {

class cantera_manager
{

 public:
  cantera_manager(const zerork::mechanism& mech);
  ~cantera_manager() ;
  void equilibrate(const double p, double& T, double* y);
 private:
  Cantera::IdealGasPhase gas;
};

} //namespace zerork

#endif
