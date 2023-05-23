#ifndef NASA_POLY_CUDA_H
#define NASA_POLY_CUDA_H

#include "nasa_poly.h"
#include <hip/hip_runtime.h>

namespace zerork {

class nasa_poly_group_cuda : public nasa_poly_group
{
 public:
  nasa_poly_group_cuda(const int inpSpc, double const * const inpCoef, 
		  double const * const inpTlow, double const * const inpThigh);

  ~nasa_poly_group_cuda();

  void getG_RT_CUDA(const double T, double * G_RT_dev, hipStream_t stream) const;
  void getG_RT_CUDA_mr(const int nReactors, const double *T_dev, double *G_RT_dev, hipStream_t stream) const;
  void getH_RT_CUDA_mr(const int nReactors, const double *T_dev, double *H_RT_dev, hipStream_t stream) const;
  void getCp_R_CUDA_mr(const int nReactors, const double *T_dev, double *Cp_R_dev, hipStream_t stream) const;

  double *thermoCoeff_dev;
};

} // namespace zerork

#endif
