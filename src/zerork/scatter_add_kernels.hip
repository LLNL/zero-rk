#include "hip/hip_runtime.h"
#include "scatter_add_kernels.h"

#if defined(__CUDA_ARCH__) &&  __CUDA_ARCH__ < 600
namespace {
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                        __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);
    return __longlong_as_double(old);
}
}
#endif

namespace zerork {

// SIMD implementation of the scatter-add (or scatter accummulate) opeartion:
// 
//     dest[destId[j]]+=src[srcId[j]]
//
// The data in the destination and source arrays are ordered in so-called 
// data order.  This means that consecutive elements relate to separate data
// values corresponding to the same instruction.
//
//

//   outer loop (j): over the number of instructions
//   inner loop (k): over the number of datum operated on by the same instr
__global__ void scatterAdd_gpu_atomic_global(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[])
{
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  if(tid < nData)
  {
      int threadSrcId  = srcId[blockIdx.y]*nData+tid;
      int threadDestId = destId[blockIdx.y]*nData+tid;

      double val = src[threadSrcId];
      atomicAdd(&dest[threadDestId],val);
  }
}

__global__ void scatterAdd_gpu_atomic_global_2op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[])
{
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  if(tid < nData)
  {
      int threadSrcId1  = srcId[2*blockIdx.y]*nData+tid;
      int threadSrcId2  = srcId[2*blockIdx.y+1]*nData+tid;
      int threadDestId = destId[2*blockIdx.y]*nData+tid;
  
      double val = src[threadSrcId1]+src[threadSrcId2];
      atomicAdd(&dest[threadDestId],val);
  }
}

__global__ void scatterAdd_gpu_atomic_global_4op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[])
{
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  if(tid < nData)
  {
      int threadSrcId1  = srcId[4*blockIdx.y]*nData+tid;
      int threadSrcId2  = srcId[4*blockIdx.y+1]*nData+tid;
      int threadSrcId3  = srcId[4*blockIdx.y+2]*nData+tid;
      int threadSrcId4  = srcId[4*blockIdx.y+3]*nData+tid;
      int threadDestId = destId[4*blockIdx.y]*nData+tid;
  
      double val = src[threadSrcId1]+
                   src[threadSrcId2]+
                   src[threadSrcId3]+
                   src[threadSrcId4];
      atomicAdd(&dest[threadDestId],val);
  }
}

__global__ void scatterAdd_gpu_atomic_global_8op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[])
{
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  if(tid < nData)
  {
      __shared__ int threadSrcId[8];
      int threadDestId;
  
      int counter = threadIdx.x;
      int stride = std::min((unsigned int) blockDim.x,nData-blockDim.x*blockIdx.x);
      while(counter < 8)
      {
          threadSrcId[counter] = srcId[8*blockIdx.y+counter]*nData;
          counter += stride;
      }
      threadDestId = destId[8*blockIdx.y]*nData+tid;

      __syncthreads();

      double val = src[threadSrcId[0]+tid]+
                   src[threadSrcId[1]+tid]+
                   src[threadSrcId[2]+tid]+
                   src[threadSrcId[3]+tid]+
                   src[threadSrcId[4]+tid]+
                   src[threadSrcId[5]+tid]+
                   src[threadSrcId[6]+tid]+
                   src[threadSrcId[7]+tid];
      atomicAdd(&dest[threadDestId],val);
  }
}

__global__ void scatterAdd_gpu_atomic_global_16op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[])
{
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  if(tid < nData)
  {
      __shared__ int threadSrcId[16];
      int threadDestId;

      int counter = threadIdx.x;
      int stride = std::min((unsigned int) blockDim.x,nData-blockDim.x*blockIdx.x);
      while(counter < 16)
      {
          threadSrcId[counter] = srcId[16*blockIdx.y+counter]*nData;
          counter += stride;
      }
      threadDestId = destId[16*blockIdx.y]*nData+tid;

      __syncthreads();

      double val = src[threadSrcId[0 ]+tid]+
                   src[threadSrcId[1 ]+tid]+
                   src[threadSrcId[2 ]+tid]+
                   src[threadSrcId[3 ]+tid]+
                   src[threadSrcId[4 ]+tid]+
                   src[threadSrcId[5 ]+tid]+
                   src[threadSrcId[6 ]+tid]+
                   src[threadSrcId[7 ]+tid]+
                   src[threadSrcId[8 ]+tid]+
                   src[threadSrcId[9 ]+tid]+
                   src[threadSrcId[10]+tid]+
                   src[threadSrcId[11]+tid]+
                   src[threadSrcId[12]+tid]+
                   src[threadSrcId[13]+tid]+
                   src[threadSrcId[14]+tid]+
                   src[threadSrcId[15]+tid];
      atomicAdd(&dest[threadDestId],val);
  }
}

__global__ void scatterAdd_gpu_atomic_global_32op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[])
{
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  if(tid < nData)
  {
      __shared__ int threadSrcId[32];
      int threadDestId;
  
      int counter = threadIdx.x;
      int stride = std::min((unsigned int) blockDim.x,nData-blockDim.x*blockIdx.x);
      while(counter < 32)
      {
          threadSrcId[counter] = srcId[32*blockIdx.y+counter]*nData;
          counter += stride;
      }
      threadDestId = destId[32*blockIdx.y]*nData+tid;

      __syncthreads();
 
      double val = src[threadSrcId[0 ]+tid]+
                   src[threadSrcId[1 ]+tid]+
                   src[threadSrcId[2 ]+tid]+
                   src[threadSrcId[3 ]+tid]+
                   src[threadSrcId[4 ]+tid]+
                   src[threadSrcId[5 ]+tid]+
                   src[threadSrcId[6 ]+tid]+
                   src[threadSrcId[7 ]+tid]+
                   src[threadSrcId[8 ]+tid]+
                   src[threadSrcId[9 ]+tid]+
                   src[threadSrcId[10]+tid]+
                   src[threadSrcId[11]+tid]+
                   src[threadSrcId[12]+tid]+
                   src[threadSrcId[13]+tid]+
                   src[threadSrcId[14]+tid]+
                   src[threadSrcId[15]+tid]+
                   src[threadSrcId[16]+tid]+
                   src[threadSrcId[17]+tid]+
                   src[threadSrcId[18]+tid]+
                   src[threadSrcId[19]+tid]+
                   src[threadSrcId[20]+tid]+
                   src[threadSrcId[21]+tid]+
                   src[threadSrcId[22]+tid]+
                   src[threadSrcId[23]+tid]+
                   src[threadSrcId[24]+tid]+
                   src[threadSrcId[25]+tid]+
                   src[threadSrcId[26]+tid]+
                   src[threadSrcId[27]+tid]+
                   src[threadSrcId[28]+tid]+
                   src[threadSrcId[29]+tid]+
                   src[threadSrcId[30]+tid]+
                   src[threadSrcId[31]+tid];
      atomicAdd(&dest[threadDestId],val);
  }
}

__global__ void scatterAdd_gpu_atomic_global_64op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[])
{
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  if(tid < nData)
  {
      __shared__ int threadSrcId[64];
      int threadDestId;

      int counter = threadIdx.x;
      int stride = std::min((unsigned int) blockDim.x,nData-blockDim.x*blockIdx.x);
      while(counter < 64)
      {
          threadSrcId[counter] = srcId[64*blockIdx.y+counter]*nData;
          counter += stride;
      }
      threadDestId = destId[64*blockIdx.y]*nData+tid;

      __syncthreads();

      double val = 0;
      #pragma unroll 32
      for (size_t idx = 0; idx < 64; ++idx)
	  val += src[threadSrcId[idx]+tid];
      atomicAdd(&dest[threadDestId],val);
  }
}

__global__ void scatterAdd_gpu_atomic_global_128op(const int nOps,
				      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[])
{
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  if(tid < nData)
  {
      __shared__ int threadSrcId[128];
      int threadDestId;

      int counter = threadIdx.x;
      int stride = std::min((unsigned int) blockDim.x,nData-blockDim.x*blockIdx.x);
      while(counter < 128)
      {
          threadSrcId[counter] = srcId[128*blockIdx.y+counter]*nData;
          counter += stride;
      }
      threadDestId = destId[128*blockIdx.y]*nData+tid;

      __syncthreads();

      double val = 0;
      #pragma unroll 32
      for (size_t idx = 0; idx < 128; ++idx)
	  val += src[threadSrcId[idx]+tid];
      atomicAdd(&dest[threadDestId],val);
  }
}

__global__ void scatterAdd_gpu_atomic_global_fused(
                                      const int nOps[],
                                      const int srcId[],
				      const int destId[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[])
{
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  if(tid < nData)
  {
    int nAdd = 128;
    __shared__ int threadSrcId[128];
    for(int j = 7; j >= 0; j--) {
      int nOpsCurr = (nOps[j] - nOps[j+1])/nAdd;
      int threadDestId;

      if(blockIdx.y < nOpsCurr) {
        const int* srcIdCurr = &srcId[nOps[j+1]];
        const int* destIdCurr = &destId[nOps[j+1]];

        int counter = threadIdx.x;
        int stride = std::min((unsigned int) blockDim.x,nData-blockDim.x*blockIdx.x);
        while(counter < nAdd)
        {
            threadSrcId[counter] = srcIdCurr[nAdd*blockIdx.y+counter]*nData;
            counter += stride;
        }
        threadDestId = destIdCurr[nAdd*blockIdx.y]*nData+tid;
      }

      __syncthreads();

      if(blockIdx.y < nOpsCurr) {
        double val = 0.0;
        for(int k = 0; k < nAdd; k++) {
          val += src[threadSrcId[k]+tid];
        }
        atomicAdd(&dest[threadDestId],val);
      }
      __syncthreads();
      nAdd /= 2;
    }
  }
}

__global__ void multScatterAdd_gpu_atomic_global_fused(
                                      const int nOps[],
                                      const int srcId[],
				      const int destId[], 
				      const double srcMult[], 
				      const int nData,
				      const int srcSize,
				      const double src[],
				      const int destSize,
				      double dest[])
{
  int tid = blockDim.x*blockIdx.x+threadIdx.x;
  if(tid < nData)
  {
    int nAdd = 128;
    __shared__ int threadSrcId[128];
    for(int j = 7; j >= 0; j--) {
      int nOpsCurr = (nOps[j] - nOps[j+1])/nAdd;
      int threadDestId;

      if(blockIdx.y < nOpsCurr) {
        const int* srcIdCurr = &srcId[nOps[j+1]];
        const int* destIdCurr = &destId[nOps[j+1]];

        int counter = threadIdx.x;
        int stride = std::min((unsigned int) blockDim.x,nData-blockDim.x*blockIdx.x);
        while(counter < nAdd)
        {
            threadSrcId[counter] = srcIdCurr[nAdd*blockIdx.y+counter];
            counter += stride;
        }
        threadDestId = destIdCurr[nAdd*blockIdx.y]*nData+tid;
      }

      __syncthreads();

      if(blockIdx.y < nOpsCurr) {
        double val = 0.0;
        for(int k = 0; k < nAdd; k++) {
          val += srcMult[nOps[j+1]+nAdd*blockIdx.y+k]*src[threadSrcId[k]*nData+tid];
        }
        atomicAdd(&dest[threadDestId],val);
      }
      __syncthreads();
      nAdd /= 2;
    }
  }
}

} // namespace zerork

