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
  __shared__ int threadSrcId[8];
  if(tid < nData)
  {
      int counter = threadIdx.x;
      int stride = min(blockDim.x,nData-blockDim.x*blockIdx.x);
      while(counter < 8)
      {
          threadSrcId[counter] = srcId[8*blockIdx.y+counter]*nData;
          counter += stride;
      }
  }

  __syncthreads();

  if(tid < nData)
  {
      int threadDestId = destId[8*blockIdx.y]*nData+tid;
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
  __shared__ int threadSrcId[16];
  if(tid < nData)
  {
      int counter = threadIdx.x;
      int stride = min(blockDim.x,nData-blockDim.x*blockIdx.x);
      while(counter < 16)
      {
          threadSrcId[counter] = srcId[16*blockIdx.y+counter]*nData;
          counter += stride;
      }
  }

  __syncthreads();

  if(tid < nData)
  {
      int threadDestId = destId[16*blockIdx.y]*nData+tid;
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
  __shared__ int threadSrcId[32];
  if(tid < nData)
  {
      int counter = threadIdx.x;
      int stride = min(blockDim.x,nData-blockDim.x*blockIdx.x);
      while(counter < 32)
      {
          threadSrcId[counter] = srcId[32*blockIdx.y+counter]*nData;
          counter += stride;
      }
  }

  __syncthreads();

  if(tid < nData)
  {
      int threadDestId = destId[32*blockIdx.y]*nData+tid;
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
  __shared__ int threadSrcId[64];
  if(tid < nData)
  {
      int counter = threadIdx.x;
      int stride = min(blockDim.x,nData-blockDim.x*blockIdx.x);
      while(counter < 64)
      {
          threadSrcId[counter] = srcId[64*blockIdx.y+counter]*nData;
          counter += stride;
      }
  }

  __syncthreads();

  if(tid < nData)
  {
      int threadDestId = destId[64*blockIdx.y]*nData+tid;
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
                   src[threadSrcId[31]+tid]+
                   src[threadSrcId[32]+tid]+
                   src[threadSrcId[33]+tid]+
                   src[threadSrcId[34]+tid]+
                   src[threadSrcId[35]+tid]+
                   src[threadSrcId[36]+tid]+
                   src[threadSrcId[37]+tid]+
                   src[threadSrcId[38]+tid]+
                   src[threadSrcId[39]+tid]+
                   src[threadSrcId[40]+tid]+
                   src[threadSrcId[41]+tid]+
                   src[threadSrcId[42]+tid]+
                   src[threadSrcId[43]+tid]+
                   src[threadSrcId[44]+tid]+
                   src[threadSrcId[45]+tid]+
                   src[threadSrcId[46]+tid]+
                   src[threadSrcId[47]+tid]+
                   src[threadSrcId[48]+tid]+
                   src[threadSrcId[49]+tid]+
                   src[threadSrcId[50]+tid]+
                   src[threadSrcId[51]+tid]+
                   src[threadSrcId[52]+tid]+
                   src[threadSrcId[53]+tid]+
                   src[threadSrcId[54]+tid]+
                   src[threadSrcId[55]+tid]+
                   src[threadSrcId[56]+tid]+
                   src[threadSrcId[57]+tid]+
                   src[threadSrcId[58]+tid]+
                   src[threadSrcId[59]+tid]+
                   src[threadSrcId[60]+tid]+
                   src[threadSrcId[61]+tid]+
                   src[threadSrcId[62]+tid]+
                   src[threadSrcId[63]+tid];
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
  __shared__ int threadSrcId[128];
  if(tid < nData)
  {
      int counter = threadIdx.x;
      int stride = min(blockDim.x,nData-blockDim.x*blockIdx.x);
      while(counter < 128)
      {
          threadSrcId[counter] = srcId[128*blockIdx.y+counter]*nData;
          counter += stride;
      }
  }

  __syncthreads();

  if(tid < nData)
  {
      int threadDestId = destId[128*blockIdx.y]*nData+tid;
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
                   src[threadSrcId[31]+tid]+
                   src[threadSrcId[32]+tid]+
                   src[threadSrcId[33]+tid]+
                   src[threadSrcId[34]+tid]+
                   src[threadSrcId[35]+tid]+
                   src[threadSrcId[36]+tid]+
                   src[threadSrcId[37]+tid]+
                   src[threadSrcId[38]+tid]+
                   src[threadSrcId[39]+tid]+
                   src[threadSrcId[40]+tid]+
                   src[threadSrcId[41]+tid]+
                   src[threadSrcId[42]+tid]+
                   src[threadSrcId[43]+tid]+
                   src[threadSrcId[44]+tid]+
                   src[threadSrcId[45]+tid]+
                   src[threadSrcId[46]+tid]+
                   src[threadSrcId[47]+tid]+
                   src[threadSrcId[48]+tid]+
                   src[threadSrcId[49]+tid]+
                   src[threadSrcId[50]+tid]+
                   src[threadSrcId[51]+tid]+
                   src[threadSrcId[52]+tid]+
                   src[threadSrcId[53]+tid]+
                   src[threadSrcId[54]+tid]+
                   src[threadSrcId[55]+tid]+
                   src[threadSrcId[56]+tid]+
                   src[threadSrcId[57]+tid]+
                   src[threadSrcId[58]+tid]+
                   src[threadSrcId[59]+tid]+
                   src[threadSrcId[60]+tid]+
                   src[threadSrcId[61]+tid]+
                   src[threadSrcId[62]+tid]+
                   src[threadSrcId[63]+tid]+
                   src[threadSrcId[64]+tid]+
                   src[threadSrcId[65]+tid]+
                   src[threadSrcId[66]+tid]+
                   src[threadSrcId[67]+tid]+
                   src[threadSrcId[68]+tid]+
                   src[threadSrcId[69]+tid]+
                   src[threadSrcId[70]+tid]+
                   src[threadSrcId[71]+tid]+
                   src[threadSrcId[72]+tid]+
                   src[threadSrcId[73]+tid]+
                   src[threadSrcId[74]+tid]+
                   src[threadSrcId[75]+tid]+
                   src[threadSrcId[76]+tid]+
                   src[threadSrcId[77]+tid]+
                   src[threadSrcId[78]+tid]+
                   src[threadSrcId[79]+tid]+
                   src[threadSrcId[80]+tid]+
                   src[threadSrcId[81]+tid]+
                   src[threadSrcId[82]+tid]+
                   src[threadSrcId[83]+tid]+
                   src[threadSrcId[84]+tid]+
                   src[threadSrcId[85]+tid]+
                   src[threadSrcId[86]+tid]+
                   src[threadSrcId[87]+tid]+
                   src[threadSrcId[88]+tid]+
                   src[threadSrcId[89]+tid]+
                   src[threadSrcId[90]+tid]+
                   src[threadSrcId[91]+tid]+
                   src[threadSrcId[92]+tid]+
                   src[threadSrcId[93]+tid]+
                   src[threadSrcId[94]+tid]+
                   src[threadSrcId[95]+tid]+
                   src[threadSrcId[96]+tid]+
                   src[threadSrcId[97]+tid]+
                   src[threadSrcId[98]+tid]+
                   src[threadSrcId[99]+tid]+
                   src[threadSrcId[100]+tid]+
                   src[threadSrcId[101]+tid]+
                   src[threadSrcId[102]+tid]+
                   src[threadSrcId[103]+tid]+
                   src[threadSrcId[104]+tid]+
                   src[threadSrcId[105]+tid]+
                   src[threadSrcId[106]+tid]+
                   src[threadSrcId[107]+tid]+
                   src[threadSrcId[108]+tid]+
                   src[threadSrcId[109]+tid]+
                   src[threadSrcId[110]+tid]+
                   src[threadSrcId[111]+tid]+
                   src[threadSrcId[112]+tid]+
                   src[threadSrcId[113]+tid]+
                   src[threadSrcId[114]+tid]+
                   src[threadSrcId[115]+tid]+
                   src[threadSrcId[116]+tid]+
                   src[threadSrcId[117]+tid]+
                   src[threadSrcId[118]+tid]+
                   src[threadSrcId[119]+tid]+
                   src[threadSrcId[120]+tid]+
                   src[threadSrcId[121]+tid]+
                   src[threadSrcId[122]+tid]+
                   src[threadSrcId[123]+tid]+
                   src[threadSrcId[124]+tid]+
                   src[threadSrcId[125]+tid]+
                   src[threadSrcId[126]+tid]+
                   src[threadSrcId[127]+tid];
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
  __shared__ int threadSrcId[128];
  int nAdd = 128;
  for(int j = 7; j >= 0; j--) {
    int nOpsCurr = (nOps[j] - nOps[j+1])/nAdd;
    const int* srcIdCurr = &srcId[nOps[j+1]];
    const int* destIdCurr = &destId[nOps[j+1]];
    if(tid < nData)
    {
      if(blockIdx.y < nOpsCurr) {
        int counter = threadIdx.x;
        int stride = min(blockDim.x,nData-blockDim.x*blockIdx.x);
        while(counter < nAdd)
        {
            threadSrcId[counter] = srcIdCurr[nAdd*blockIdx.y+counter]*nData;
            counter += stride;
        }
      }
    }

    __syncthreads();

    if(tid < nData)
    {
      if(blockIdx.y < nOpsCurr) {
        int threadDestId = destIdCurr[nAdd*blockIdx.y]*nData+tid;
        double val = 0.0;
        for(int k = 0; k < nAdd; k++) {
          val += src[threadSrcId[k]+tid];
        }
        atomicAdd(&dest[threadDestId],val);
      }
    }
    __syncthreads();
    nAdd /= 2;
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
  __shared__ int threadSrcId[128];
  int nAdd = 128;
  for(int j = 7; j >= 0; j--) {
    int nOpsCurr = (nOps[j] - nOps[j+1])/nAdd;
    const int* srcIdCurr = &srcId[nOps[j+1]];
    const int* destIdCurr = &destId[nOps[j+1]];
    if(tid < nData)
    {
      if(blockIdx.y < nOpsCurr) {
        int counter = threadIdx.x;
        int stride = min(blockDim.x,nData-blockDim.x*blockIdx.x);
        while(counter < nAdd)
        {
            threadSrcId[counter] = srcIdCurr[nAdd*blockIdx.y+counter];
            counter += stride;
        }
      }
    }

    __syncthreads();

    if(tid < nData)
    {
      if(blockIdx.y < nOpsCurr) {
        int threadDestId = destIdCurr[nAdd*blockIdx.y]*nData+tid;
        double val = 0.0;
        for(int k = 0; k < nAdd; k++) {
          val += srcMult[nOps[j+1]+nAdd*blockIdx.y+k]*src[threadSrcId[k]*nData+tid];
        }
        atomicAdd(&dest[threadDestId],val);
      }
    }
    __syncthreads();
    nAdd /= 2;
  }
}

} // namespace zerork

