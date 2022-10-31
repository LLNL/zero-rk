#ifndef REACTOR_CONSTANT_PRESSURE_CPU_H
#define REACTOR_CONSTANT_PRESSURE_CPU_H

#include "reactor_nvector_serial.h"

class ReactorConstantPressureCPU : public ReactorNVectorSerial
{
 public:
  ReactorConstantPressureCPU(std::shared_ptr<zerork::mechanism> mech_ptr);
  ~ReactorConstantPressureCPU();

  void InitializeState(const double reactor_time,
                       const int n_reactors,
                       const double *T,
                       const double *P,
                       const double *mf,
                       const double *dpdt,
                       const double *e_src,
                       const double *y_src);

  void GetState(const double reactor_time,
                double *T,
                double *P,
                double *mf);

  int GetTimeDerivative(const double reactor_time,
                        N_Vector state,
                        N_Vector derivative);

};

#endif
