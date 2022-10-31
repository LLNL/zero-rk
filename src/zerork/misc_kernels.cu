#include "misc_kernels.h"
#include "zerork_cuda_defs.h"
#include "constants.h"
#include <stdio.h>
#include <algorithm> //std::min


namespace zerork {

// Kernels:

void __global__ kernel_printCudaArray(const int n, const double *A_dev)
{
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  int stride = gridDim.x*blockDim.x;
  while(tid < n)
  {
      printf("%d : %g\n",tid,A_dev[tid]);
      tid+=stride;
  }
}

void __global__ kernel_setDoubleArrayVal(double *array, const double  val, const int idxStart, const int idxEnd)
{
    int tidx = threadIdx.x + blockDim.x * blockIdx.x + idxStart;
    int stride = blockDim.x * gridDim.x;

    for(; tidx < idxEnd; tidx += stride)
        array[tidx] = val;
}

// Wrappers:

void printCudaArray(const int n, const double *A_dev)
{
  int nThreads = MAX_THREADS_PER_BLOCK;
  int nBlocks = (n + nThreads -1 )/nThreads;
  kernel_printCudaArray<<<nBlocks, nThreads>>>(n,A_dev);
#ifdef ZERORK_FULL_DEBUG
  cudaDeviceSynchronize();
  checkCudaError(cudaGetLastError(),"printCudaArray()");
#endif
}

void setDoubleArrayVal(double * array,const double value,const int idxStart,const int idxEnd, const cudaStream_t stream)
{
  int nVals = idxEnd - idxStart; //not inclusive
  int nThreads = std::min(nVals,MAX_THREADS_PER_BLOCK);
  int nBlocks = (nVals+nThreads-1)/nThreads;
  kernel_setDoubleArrayVal<<<nBlocks,nThreads,0,stream>>>(array,value,idxStart,idxEnd);
#ifdef ZERORK_FULL_DEBUG
  cudaDeviceSynchronize();
  checkCudaError(cudaGetLastError(),"setDoubleArrayVal()");
#endif

}
} // namespace zerork

