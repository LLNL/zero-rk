#ifndef REACTOR_CONSTANT_VOLUME_H
#define REACTOR_CONSTANT_VOLUME_H
#include <memory>

#include "reactor_nvector_serial.h"

class ReactorConstantVolume : public ReactorNVectorSerial
{
 public:
  ReactorConstantVolume(std::shared_ptr<zerork::mechanism> mech_ptr);
  ~ReactorConstantVolume();

  void InitializeState(const double reactor_time,
                       const int n_reactors,
                       const double *T,
                       const double *P,
                       const double *mf,
                       const double *dpdt,
                       const double *e_src,
                       const double *y_src);

  void GetState(double *T,
                double *P,
                double *mf);

  int GetTimeDerivative(const double reactor_time,
                        N_Vector state,
                        N_Vector derivative);

};

#endif
