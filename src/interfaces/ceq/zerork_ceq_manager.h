
#ifndef ZERORK_CEQ_MANAGER_H
#define ZERORK_CEQ_MANAGER_H

#include "mechanism.h"

namespace zerork {

class ceq_manager
{

 public:
  ceq_manager(const zerork::mechanism& mech);
  ~ceq_manager() ;
  void equilibrate(const double p, double& T, double* y);
};

} //namespace zerork

#endif
