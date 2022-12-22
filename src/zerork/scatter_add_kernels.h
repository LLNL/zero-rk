#ifndef ZERORK_SCATTER_ADD_KERNELS_H
#define ZERORK_SCATTER_ADD_KERNELS_H

namespace zerork {

__global__ void scatterAdd_gpu_atomic_global(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[]);
__global__ void scatterAdd_gpu_atomic_global_2op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[]);
__global__ void scatterAdd_gpu_atomic_global_4op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[]);
__global__ void scatterAdd_gpu_atomic_global_8op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[]);
__global__ void scatterAdd_gpu_atomic_global_16op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[]);
__global__ void scatterAdd_gpu_atomic_global_32op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[]);
__global__ void scatterAdd_gpu_atomic_global_64op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[]);
__global__ void scatterAdd_gpu_atomic_global_128op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[]);
__global__ void scatterAdd_gpu_atomic_global_fused(
                                      const int nOps[],
                                      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[]);

__global__ void multScatterAdd_gpu_atomic_global_fused(
                                      const int nOps[],
                                      const int srcId[],
				      const int destId[], 
				      const double srcMult[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[]);

} // namespace zerork

#endif
