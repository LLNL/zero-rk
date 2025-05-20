#ifndef NASA_POLY_KERNELS_H
#define NASA_POLY_KERNELS_H

namespace zerork {

void nasa_poly_group_getG_RT_CUDA_mr(const int nReactors, const double *T_dev, double *G_RT_dev, const int nGroupSpc, const double * thermoCoeff_dev, cudaStream_t stream);
void nasa_poly_group_getH_RT_CUDA_mr(const int nReactors, const double *T_dev, double *H_RT_dev, const int nGroupSpc, const double * thermoCoeff_dev, cudaStream_t stream);
void nasa_poly_group_getCp_R_CUDA_mr(const int nReactors, const double *T_dev, double *Cp_R_dev, const int nGroupSpc, const double * thermoCoeff_dev, cudaStream_t stream);

} // namespace zerork

#endif
