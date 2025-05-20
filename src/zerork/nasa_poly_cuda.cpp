#include "nasa_poly_cuda.h"
#include "constants.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <cuda_runtime.h>
#include "zerork_cuda_defs.h"
#include "nasa_poly_kernels.h"

namespace zerork {

nasa_poly_group_cuda::nasa_poly_group_cuda(const int inpSpc,
  double const * const inpCoef, double const * const inpTlow,
  double const * const inpThigh)
    : nasa_poly_group(inpSpc,inpCoef,inpTlow,inpThigh)
{
  checkCudaError
  (
     cudaMalloc((void**)&thermoCoeff_dev,sizeof(double)*inpSpc*LDA_THERMO_POLY_D5R2),
     "cudaMalloc(... thermoCoeff_dev ...)"
  );
  cudaMemcpy(thermoCoeff_dev,thermoCoef,inpSpc*LDA_THERMO_POLY_D5R2*sizeof(double),cudaMemcpyHostToDevice);

#ifdef ZERORK_FULL_DEBUG
  cudaDeviceSynchronize();
  checkCudaError(cudaGetLastError(),"nasa_poly_group::nasa_poly_group( ... )");
#endif
}

nasa_poly_group_cuda::~nasa_poly_group_cuda()
{
  cudaFree(thermoCoeff_dev);
}

void nasa_poly_group_cuda::getG_RT_CUDA_mr(const int nReactors, const double *T_dev, double *G_RT_dev, cudaStream_t stream) const
{
    nasa_poly_group_getG_RT_CUDA_mr(nReactors, T_dev, G_RT_dev, nGroupSpc, thermoCoeff_dev, stream);
}

void nasa_poly_group_cuda::getH_RT_CUDA_mr(const int nReactors, const double *T_dev, double *H_RT_dev, cudaStream_t stream) const
{
    nasa_poly_group_getH_RT_CUDA_mr(nReactors, T_dev, H_RT_dev, nGroupSpc, thermoCoeff_dev, stream);
}

void nasa_poly_group_cuda::getCp_R_CUDA_mr(const int nReactors, const double *T_dev, double *Cp_R_dev, cudaStream_t stream) const
{
    nasa_poly_group_getCp_R_CUDA_mr(nReactors, T_dev, Cp_R_dev, nGroupSpc, thermoCoeff_dev, stream);
}

} // namespace zerork

