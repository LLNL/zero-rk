#include <algorithm> //std::min
#include "mechanism_kernels.h"
#include "zerork_cuda_defs.h"
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

void __global__ kernel_minusEqualOne(const int N, double * A_dev)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  if(tid < N)
  {
      A_dev[tid] -= 1.0;
  }
}

void __global__ kernel_invert(const int N, double * A_dev)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  if(tid < N)
  {
      A_dev[tid] = 1.0/A_dev[tid];
  }
}

void __global__ kernel_meanCpMass_mr
(
    const int nReactors,
    const int nSpc,
    const double *RuInvMolWt_dev,
    const double *y_dev,
    double *Cp_dev, // species Cp_dev
    double *meanCp_dev // mean Cp for reactor
)
{
    int reactorid = blockIdx.x*blockDim.x + threadIdx.x;
    int tid = blockIdx.y*blockDim.y + threadIdx.y;
    extern __shared__ double values[];
    if(tid < nSpc)
    {
        if(reactorid < nReactors)
        {
            int valIdx = blockDim.x*threadIdx.y+threadIdx.x;
            values[valIdx] = Cp_dev[nReactors*tid+reactorid];
            values[valIdx] = values[valIdx]*RuInvMolWt_dev[tid];
            Cp_dev[nReactors*tid+reactorid] = values[valIdx];
            values[valIdx] *= y_dev[nReactors*tid+reactorid];
	}
    }
    __syncthreads();
    if(tid < nSpc)
    {
        if(reactorid < nReactors)
        {
            if(threadIdx.y == 0)
            {
                double accum = 0.0;
                int lim = min(blockDim.y,nSpc-blockDim.y*blockIdx.y);
                for(int i = 0; i < lim; ++i)
                {
                    accum += values[blockDim.x*i+threadIdx.x];
                }
                atomicAdd(&meanCp_dev[reactorid],accum);
            }
        }
    }
}

void __global__ kernel_Cv_from_Cp_mr
(
    const int nReactors,
    const int nSpc,
    const double *RuInvMolWt_dev,
    const double *y_dev,
    double *Cv_dev, //On entry Cp_dev, on exit Cv_dev
    double *meanCv_dev // mean Cv for reactor
)
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    int stride = gridDim.x*blockDim.x;
    while(tid < nSpc*nReactors)
    {
       int speciesid = tid / nReactors;
       int reactorid = tid % nReactors;
       double val = Cv_dev[tid];
       val = (val-1.0)*RuInvMolWt_dev[speciesid];
       Cv_dev[tid] = val;
       val *= y_dev[tid];
       atomicAdd(&meanCv_dev[reactorid],val);
       tid += stride;
    }
}

void __global__ kernel_getCfromVY_mr_dev
(
    const int nReactors,
    const int nSpc,
    const double *v_dev,
    const double *y_dev,
    const double *invMolWt_dev,
    double *c_dev
)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = gridDim.x*blockDim.x;
  int speciesid = tid / nReactors;
  int reactorid = tid % nReactors;
  while(tid < nSpc*nReactors)
  {
    double dens = 1.0/v_dev[reactorid];
    c_dev[tid] = dens*invMolWt_dev[speciesid]*y_dev[tid];
    tid += stride;
  }
}

void __global__ kernel_getMassIntEnergyFromTY_mr_dev
(
    const int nReactors,
    const int nSpc,
    const double *T_dev,
    const double *y_dev,
    const double *RuInvMolWt_dev,
    double *u_spc_dev,
    double *u_dev
)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int reactorid = blockIdx.y*blockDim.y + threadIdx.y;
  extern __shared__ double values[];
  if(tid < nSpc)
  {
      if(reactorid < nReactors)
      {
            int valIdx = blockDim.y*threadIdx.x+threadIdx.y;
            u_spc_dev[nReactors*tid+reactorid]  = (u_spc_dev[nReactors*tid+reactorid] - 1)*RuInvMolWt_dev[tid]*T_dev[reactorid];
            values[valIdx] = u_spc_dev[nReactors*tid+reactorid]*y_dev[nReactors*tid+reactorid];
      }
  }
  __syncthreads();
  if(tid < nSpc)
  {
      if(reactorid < nReactors)
      {
            if(threadIdx.x == 0)
            {
                double accum = 0.0;
                int lim = min(blockDim.x,nSpc-blockDim.x*blockIdx.x);
                for(int i = 0; i < lim; ++i)
                {
                    accum += values[blockDim.y*i+threadIdx.y];
                }
                atomicAdd(&u_dev[reactorid],accum);
            }
      }
  }
}

void __global__ kernel_getMassEnthalpyFromTY_mr_dev
(
    const int nReactors,
    const int nSpc,
    const double *T_dev,
    const double *y_dev,
    const double *RuInvMolWt_dev,
    double *h_spc_dev,
    double *h_dev
)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int reactorid = blockIdx.y*blockDim.y + threadIdx.y;
  extern __shared__ double values[];
  if(tid < nSpc)
  {
      if(reactorid < nReactors)
      {
            int valIdx = blockDim.y*threadIdx.x+threadIdx.y;
            h_spc_dev[nReactors*tid+reactorid]  = h_spc_dev[nReactors*tid+reactorid]*RuInvMolWt_dev[tid]*T_dev[reactorid];
            values[valIdx] = h_spc_dev[nReactors*tid+reactorid]*y_dev[nReactors*tid+reactorid];
      }
  }
  __syncthreads();
  if(tid < nSpc)
  {
      if(reactorid < nReactors)
      {
            if(threadIdx.x == 0)
            {
                double accum = 0.0;
                int lim = min(blockDim.x,nSpc-blockDim.x*blockIdx.x);
                for(int i = 0; i < lim; ++i)
                {
                    accum += values[blockDim.y*i+threadIdx.y];
                }
                atomicAdd(&h_dev[reactorid],accum);
            }
      }
  }
}

void __global__ kernel_getMolWtMixFromY
(
    const int nReactors,
    const int nSpc,
    const double *y_dev,
    const double *invMolWt_dev,
    double *molWtMix_dev
)
{
  int reactorid = blockIdx.x*blockDim.x + threadIdx.x;
  int tid = blockIdx.y*blockDim.y + threadIdx.y;
  extern __shared__ double values[];
  if(tid < nSpc)
  {
      if(reactorid < nReactors)
      {
            int valIdx = blockDim.x*threadIdx.y+threadIdx.x;
            values[valIdx] = y_dev[nReactors*tid+reactorid];
            values[valIdx] = values[valIdx]*invMolWt_dev[tid];
      }
  }
  __syncthreads();
  if(tid < nSpc)
  {
      if(reactorid < nReactors)
      {
            if(threadIdx.y == 0)
            {
                double accum = 0.0;
                int lim = min(blockDim.y,nSpc-blockDim.y*blockIdx.y);
                for(int i = 0; i < lim; ++i)
                {
                    accum += values[blockDim.x*i+threadIdx.x];
                }
                atomicAdd(&molWtMix_dev[reactorid],accum);
            }
      }
  }
}


void __global__ kernel_getDensityFromTPMW(const int nReactors, const double Ru,
                                          const double T_dev[],
                                          const double P_dev[],
                                          double  dens_dev[])
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  if(tid < nReactors)
  {
      dens_dev[tid] = P_dev[tid]*dens_dev[tid]/(Ru*T_dev[tid]);
  }
}


void __global__  kernel_getTemperatureFromEY_mr_dev_part1(const int nReactors,
                                               const double T_min, const double T_max,
                                               const double* E,
                                               const double* Emin, const double* Emax,
                                               const double* cv_min, const double* cv_max,
                                               double* T, int* converged)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = gridDim.x*blockDim.x;
  while(tid < nReactors) {
    double E_local = E[tid];
    double Emin_local = Emin[tid];
    double Emax_local = Emax[tid];
    int converged_local = 0;
    if(E_local < Emin_local) {
      T[tid] = T_min - (Emin_local-E_local)/cv_min[tid];
      converged_local = 1;
    }
    if(E_local > Emax_local) {
      T[tid] = T_max - (Emax_local-E_local)/cv_max[tid];
      converged_local = 1;
    }
    if(converged_local == 0) {
      double T_local = T[tid];
      if(T_local < T_min || T_local > T_max) {
        T[tid]= T_min + (T_max-T_min)/(Emax_local - Emin_local)*(E_local-Emin_local);
      }
    }
    converged[tid] = converged_local;
    tid += stride;
  }
}
     

void __global__ kernel_getTemperatureFromEY_mr_dev_iter(const int nReactors, const double tolerance,
                                               const double* E, const double* Eiter, const double* cv_iter,
                                               double* T, int* converged)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = gridDim.x*blockDim.x;
  while(tid < nReactors) {
    int converged_local = converged[tid];
    if(converged_local == 0) {
      double E_local = E[tid];
      double Eiter_local = Eiter[tid];
      double cv_iter_local = cv_iter[tid];
      double delta_temp = (E_local-Eiter_local)/cv_iter_local;
      if(delta_temp > 100) delta_temp = 100;
      if(delta_temp < -100) delta_temp = -100;
      if(std::abs(delta_temp) < tolerance ) {
        converged[tid] = 1;
      } else {
        double T_local = T[tid];
        double T_update = T_local + delta_temp;
        if(T_update == T_local) {
          converged[tid] = 1;
        } else {
          T[tid] = T_update;
        }
      }
    }
    tid += stride;
  }
}

void __global__ kernel_getPressureFromTVY(const int nReactors, const double Ru,
                                          const double* T,
                                          const double* v,
                                          double* P) 
{
  //on entry P has average molecular weight
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = gridDim.x*blockDim.x;
  while(tid < nReactors) {
    P[tid] = Ru*T[tid]/(v[tid]*P[tid]);
    tid += stride;
  }
}

// Wrappers:

void minusEqualOne(const int N, double *A_dev)
{
  int nThreads = MAX_THREADS_PER_BLOCK;
  int nBlocks = (N + nThreads - 1)/nThreads;
  kernel_minusEqualOne<<<nBlocks,nThreads>>>(N,A_dev);
}

void invert(const int N, double *A_dev)
{
  int nThreads = MAX_THREADS_PER_BLOCK;
  int nBlocks = (N + nThreads - 1)/nThreads;
  kernel_invert<<<nBlocks,nThreads>>>(N,A_dev);
}

void meanCpMass_mr(const int nReactors, const int nSpc, const double *RuInvMolWt_dev,
    const double *y_dev, double *cpSpc_dev, double *cpReactors_dev)
{
  dim3 nThreads2D,nBlocks2D;
  nThreads2D.x = 512; //FIXME : Magic number
  nThreads2D.y = std::min(nSpc,MAX_THREADS_PER_BLOCK/((int)nThreads2D.x));
  nBlocks2D.x = (nReactors+nThreads2D.x-1)/nThreads2D.x;
  nBlocks2D.y = (nSpc+nThreads2D.y-1)/nThreads2D.y;
  size_t shrMemSize = sizeof(double)*nThreads2D.x*nThreads2D.y;
  kernel_meanCpMass_mr<<<nBlocks2D,nThreads2D,shrMemSize>>>(nReactors,nSpc,RuInvMolWt_dev,y_dev,cpSpc_dev,cpReactors_dev);
}

void Cv_from_Cp_mr(const int nReactors, const int nSpc, const double *RuInvMolWt_dev,
    const double *y_dev, double *cvSpc_dev, double *cvReactors_dev)
{
  int nThreads = MAX_THREADS_PER_BLOCK;
  int nBlocks = (nSpc*nReactors+nThreads-1)/nThreads;
  kernel_Cv_from_Cp_mr<<<nBlocks,nThreads>>>(nReactors,nSpc,RuInvMolWt_dev,y_dev,cvSpc_dev,cvReactors_dev);
}

void getCfromVY_mr_dev_wrapper(const int nReactors, const int nSpc,
                       const double *v_dev, const double *y_dev,
                       const double *invMolWt_dev, double *c_dev)
{
  int nThreads = MAX_THREADS_PER_BLOCK;
  int nBlocks = (nReactors*nSpc+nThreads-1)/nThreads;
  kernel_getCfromVY_mr_dev<<<nBlocks,nThreads>>>
  (
      nReactors,
      nSpc,
      v_dev,
      y_dev,
      invMolWt_dev,
      c_dev
  );
}

void getMassIntEnergyFromTY_mr_dev_wrapper(const int nReactors, const int nSpc,
                       const double *T_dev, const double *y_dev,
                       const double *RuInvMolWt_dev, double* u_spc_dev, double *u_dev)
{
  dim3 nBlocks2D, nThreads2D;
  nThreads2D.x = 512; //FIXME : Magic number
  nThreads2D.y = std::min(nReactors,MAX_THREADS_PER_BLOCK/((int)nThreads2D.x));
  nBlocks2D.x = (nSpc+nThreads2D.x-1)/nThreads2D.x;
  nBlocks2D.y = (nReactors+nThreads2D.y-1)/nThreads2D.y;
  size_t shrMemSize = sizeof(double)*nThreads2D.x*nThreads2D.y;
  kernel_getMassIntEnergyFromTY_mr_dev<<<nBlocks2D,nThreads2D,shrMemSize>>>
  (
      nReactors,
      nSpc,
      T_dev,
      y_dev,
      RuInvMolWt_dev,
      u_spc_dev,
      u_dev
  );
}

void getMassEnthalpyFromTY_mr_dev_wrapper(const int nReactors, const int nSpc,
                       const double *T_dev, const double *y_dev,
                       const double *RuInvMolWt_dev, double* h_spc_dev, double *h_dev)
{
  dim3 nBlocks2D, nThreads2D;
  nThreads2D.x = 512; //FIXME : Magic number
  nThreads2D.y = std::min(nReactors,MAX_THREADS_PER_BLOCK/((int)nThreads2D.x));
  nBlocks2D.x = (nSpc+nThreads2D.x-1)/nThreads2D.x;
  nBlocks2D.y = (nReactors+nThreads2D.y-1)/nThreads2D.y;
  size_t shrMemSize = sizeof(double)*nThreads2D.x*nThreads2D.y;
  kernel_getMassEnthalpyFromTY_mr_dev<<<nBlocks2D,nThreads2D,shrMemSize>>>
  (
      nReactors,
      nSpc,
      T_dev,
      y_dev,
      RuInvMolWt_dev,
      h_spc_dev,
      h_dev
  );
}

void getMolWtMixFromY_mr_dev_wrapper(const int nReactors, const int nSpc,
                       const double *y_dev, const double *invMolWt_dev,
                       double *mwMix_dev)
{
  dim3 nBlocks2D, nThreads2D;
  nThreads2D.x = 512; //FIXME : Magic number
  nThreads2D.y = std::min(nSpc,MAX_THREADS_PER_BLOCK/((int)nThreads2D.x));
  nBlocks2D.x = (nReactors+nThreads2D.x-1)/nThreads2D.x;
  nBlocks2D.y = (nSpc+nThreads2D.y-1)/nThreads2D.y;
  size_t shrMemSize = sizeof(double)*nThreads2D.x*nThreads2D.y;
  kernel_getMolWtMixFromY<<<nBlocks2D,nThreads2D,shrMemSize>>>
  (
      nReactors,
      nSpc,
      y_dev,
      invMolWt_dev,
      mwMix_dev
  );
  invert(nReactors,mwMix_dev);
}

void getDensityfromTPMW_wrapper(const int nReactors, const double Ru,
                                const double T_dev[],
                                const double P_dev[], double dens_dev[])
{
  int nThreads = MAX_THREADS_PER_BLOCK;
  int nBlocks = (nReactors + nThreads - 1)/nThreads;
  kernel_getDensityFromTPMW<<<nBlocks,nThreads>>>(nReactors,Ru,T_dev,P_dev,dens_dev);
}


void getTemperatureFromEY_mr_dev_part1_wrapper(const int nReactors,
                                               const double T_min, const double T_max,
                                               const double* E,
                                               const double* Emin, const double* Emax,
                                               const double* cv_min, const double* cv_max,
                                               double* T, int* converged)
{
  int nThreads = MAX_THREADS_PER_BLOCK;
  int nBlocks = (nReactors + nThreads - 1)/nThreads;
  kernel_getTemperatureFromEY_mr_dev_part1<<<nBlocks,nThreads>>>(nReactors,
      T_min, T_max, E, Emin, Emax, cv_min, cv_max, T, converged);
}

void getTemperatureFromEY_mr_dev_iter_wrapper(const int nReactors, double tolerance,
                                              const double* E, const double* Eiter,
                                              const double* cv_iter,
                                              double* T, int* converged)
{
  int nThreads = MAX_THREADS_PER_BLOCK;
  int nBlocks = (nReactors + nThreads - 1)/nThreads;
  kernel_getTemperatureFromEY_mr_dev_iter<<<nBlocks,nThreads>>>(nReactors, tolerance,
      E, Eiter, cv_iter, T, converged);
}

void getPressureFromTVY_wrapper(const int nReactors, const double Ru,  const double* T, const double* v, double* P)
{
  //on entry P is average molecular weight
  int nThreads = MAX_THREADS_PER_BLOCK;
  int nBlocks = (nReactors + nThreads - 1)/nThreads;
  kernel_getPressureFromTVY<<<nBlocks,nThreads>>>(nReactors, Ru, T, v, P);
}

} // namespace zerork

