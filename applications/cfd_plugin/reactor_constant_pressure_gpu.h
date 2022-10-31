#ifndef REACTOR_CONSTANT_PRESSURE_GPU_H
#define REACTOR_CONSTANT_PRESSURE_GPU_H

#include "zerork/mechanism_cuda.h"
#include "reactor_nvector_serial_cuda.h"

class ReactorConstantPressureGPU : public ReactorNVectorSerialCuda
{
 public:
  ReactorConstantPressureGPU(std::shared_ptr<zerork::mechanism_cuda> mech_ptr);
  ~ReactorConstantPressureGPU();

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

  void GetAbsoluteToleranceCorrection(N_Vector correction);

  int RootFunction(double t, N_Vector y, double *root_function);

  int GetNumRootFunctions();
};

#endif
