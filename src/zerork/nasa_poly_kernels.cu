#include <algorithm> // std::min
#include "nasa_poly_kernels.h"
#include "zerork_cuda_defs.h"
#include "constants.h"

namespace zerork {


void __global__ cuda_get_G_RT(int nSpc, const double T, const double *thermoCoeff_dev, double *G_RT_dev)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;

  if(tid < nSpc)
  {
      double Tmid,g1;
      double invT=1.0/T;
      double gMult[5];

      gMult[0]=1.0-log(T);
      gMult[1]=T;
      gMult[2]=T*gMult[1];
      gMult[3]=T*gMult[2];
      gMult[4]=T*gMult[3];

                                          // g0 = 1 - ln(T)
      gMult[1]*=-0.50000000000000000000;  // g1 = - T/2
      gMult[2]*=-0.16666666666666666667;  // g2 = - T^2/6
      gMult[3]*=-0.08333333333333333333;  // g3 = - T^3/12
      gMult[4]*=-0.05000000000000000000;  // g4 = - T^4/20

      int coefAddr=LDA_THERMO_POLY_D5R2*tid;
      Tmid=thermoCoeff_dev[coefAddr];
      if(T < Tmid)
      {
          g1=thermoCoeff_dev[coefAddr+1]*gMult[0]+
             thermoCoeff_dev[coefAddr+2]*gMult[1]+
             thermoCoeff_dev[coefAddr+3]*gMult[2]+
             thermoCoeff_dev[coefAddr+4]*gMult[3]+
             thermoCoeff_dev[coefAddr+5]*gMult[4]+
             thermoCoeff_dev[coefAddr+6]*invT-
             thermoCoeff_dev[coefAddr+7];
      }
      else
      {
          g1=thermoCoeff_dev[coefAddr+8 ]*gMult[0]+
             thermoCoeff_dev[coefAddr+9 ]*gMult[1]+
             thermoCoeff_dev[coefAddr+10]*gMult[2]+
             thermoCoeff_dev[coefAddr+11]*gMult[3]+
             thermoCoeff_dev[coefAddr+12]*gMult[4]+
             thermoCoeff_dev[coefAddr+13]*invT-
             thermoCoeff_dev[coefAddr+14];
      }

      G_RT_dev[tid] = g1;
  }
}

__device__ const double gMult[16] = {
 //Tmid, g1, g2,  g3,                     g4,                    g5    g6 g7, 
      1., 1.,-0.5,-0.16666666666666666667,-0.08333333333333333333,-0.05,1.,-1,
          1.,-0.5,-0.16666666666666666667,-0.08333333333333333333,-0.05,1.,-1,0.,
      };

void __global__ cuda_get_G_RT_v2(int nSpc, const double T, const double *thermoCoeff_dev, double *G_RT_dev)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int speciesid = tid/16;
  int myOffset = tid-16*speciesid;

  if(speciesid < nSpc)
  {
      extern double __shared__ gOut[];
      gOut[threadIdx.x] = thermoCoeff_dev[tid]*gMult[myOffset];

      __syncthreads();
      if(myOffset == 0)
      {
          if(T < gOut[threadIdx.x])
          {
              G_RT_dev[speciesid] = (((gOut[threadIdx.x+5]*T+gOut[threadIdx.x+4])*T+gOut[threadIdx.x+3])*T+gOut[threadIdx.x+2])*T
                                     + gOut[threadIdx.x+1]*(1.0-log(T))           + gOut[threadIdx.x+6]/T +gOut[threadIdx.x+7]; 
          }       
          else
          {
              G_RT_dev[speciesid] = (((gOut[threadIdx.x+12]*T+gOut[threadIdx.x+11])*T+gOut[threadIdx.x+10])*T+gOut[threadIdx.x+9])*T
                                     + gOut[threadIdx.x+8]*(1.0-log(T))           + gOut[threadIdx.x+13]/T +gOut[threadIdx.x+14]; 
          }       
      }       
  }
}

void __global__ cuda_get_G_RT_mr
(
    int nReactors, int nSpc, const double *T_dev, const double *thermoCoeff_dev, double *G_RT_dev
)
{
  int reactorid = blockIdx.x*blockDim.x + threadIdx.x;
  int speciesid = blockIdx.y;
  if(speciesid < nSpc)
  {
      if(reactorid < nReactors)
      {
          int coefAddr = LDA_THERMO_POLY_D5R2*speciesid;
          __shared__ double coeffShr[16];

          int counter = threadIdx.x;
          int stride = min(blockDim.x,nReactors-blockDim.x*blockIdx.x);
          while(counter < 16)
          {
              coeffShr[counter]=thermoCoeff_dev[coefAddr+counter];
              counter += stride;
          }
          __syncthreads();

          double Tmid,g1,g2;
          double T = T_dev[reactorid];
          double invT=1.0/T;
          double gMult[5];

          gMult[0]=1.0-log(T);
          gMult[1]=T;
          gMult[2]=T*gMult[1];
          gMult[3]=T*gMult[2];
          gMult[4]=T*gMult[3];

                                              // g0 = 1 - ln(T)
          gMult[1]*=-0.50000000000000000000;  // g1 = - T/2
          gMult[2]*=-0.16666666666666666667;  // g2 = - T^2/6
          gMult[3]*=-0.08333333333333333333;  // g3 = - T^3/12
          gMult[4]*=-0.05000000000000000000;  // g4 = - T^4/20

          Tmid=coeffShr[0];
          if(T < Tmid)
          {
          g2=coeffShr[1]*gMult[0]+
             coeffShr[2]*gMult[1]+
             coeffShr[3]*gMult[2]+
             coeffShr[4]*gMult[3]+
             coeffShr[5]*gMult[4]+
             coeffShr[6]*invT-
             coeffShr[7];

          }
          else
          {
          g2=coeffShr[8 ]*gMult[0]+
             coeffShr[9 ]*gMult[1]+
             coeffShr[10]*gMult[2]+
             coeffShr[11]*gMult[3]+
             coeffShr[12]*gMult[4]+
             coeffShr[13]*invT-
             coeffShr[14];
          }

//          if(T < Tmid)
//          {
//              g2 = g1;
//          }
          G_RT_dev[reactorid+speciesid*nReactors] = g2;
      }
  }
}

void __global__ cuda_get_H_RT_mr
(
    int nReactors, int nSpc, const double *T_dev, const double *thermoCoeff_dev, double *H_RT_dev
)
{
  int reactorid = blockIdx.x*blockDim.x + threadIdx.x;
  int speciesid = blockIdx.y;
  if(speciesid < nSpc)
  {
      if(reactorid < nReactors)
      {
          int coefAddr = LDA_THERMO_POLY_D5R2*speciesid;
          __shared__ double coeffShr[16];

          int counter = threadIdx.x;
          int stride = min(blockDim.x,nReactors-blockDim.x*blockIdx.x);
          while(counter < 16)
          {
              coeffShr[counter]=thermoCoeff_dev[coefAddr+counter];
              counter += stride;
          }
          __syncthreads();

          double Tmid,h1,h2;
          double T = T_dev[reactorid];
          double invT=1.0/T;
          double hMult[4];

          hMult[0]=T;
          hMult[1]=T*hMult[0];
          hMult[2]=T*hMult[1];
          hMult[3]=T*hMult[2];

          hMult[0]*=0.50000000000000000000; // h0 = T/2
          hMult[1]*=0.33333333333333333333; // h1 = T^2/3
          hMult[2]*=0.25000000000000000000; // h2 = T^3/4
          hMult[3]*=0.20000000000000000000; // h3 = T^4/5

          Tmid=coeffShr[0];
          h1=coeffShr[1]+
             coeffShr[2]*hMult[0]+
             coeffShr[3]*hMult[1]+
             coeffShr[4]*hMult[2]+
             coeffShr[5]*hMult[3]+
             coeffShr[6]*invT;

          h2=coeffShr[8 ]+
             coeffShr[9 ]*hMult[0]+
             coeffShr[10]*hMult[1]+
             coeffShr[11]*hMult[2]+
             coeffShr[12]*hMult[3]+
             coeffShr[13]*invT;

          if(T < Tmid)
          {
              h2 = h1;
          }
          H_RT_dev[reactorid+speciesid*nReactors] = h2;
      }
  }
}

void __global__ cuda_get_Cp_R_mr
(
    int nReactors, int nSpc, const double *T_dev, const double *thermoCoeff_dev, double *Cp_R_dev
)
{
  int reactorid = blockIdx.x*blockDim.x + threadIdx.x;
  int speciesid = blockIdx.y;
  if(speciesid < nSpc)
  {
      if(reactorid < nReactors)
      {
          int coefAddr = LDA_THERMO_POLY_D5R2*speciesid;
          __shared__ double coeffShr[16];

          int counter = threadIdx.x;
          int stride = min(blockDim.x,nReactors-blockDim.x*blockIdx.x);
          while(counter < 16)
          {
              coeffShr[counter]=thermoCoeff_dev[coefAddr+counter];
              counter += stride;
          }
          __syncthreads();

          double Tmid,cp1,cp2;
          double T = T_dev[reactorid];
          double invT=1.0/T;

          Tmid=coeffShr[0];
          cp1=   coeffShr[1]+
              T*(coeffShr[2]+
              T*(coeffShr[3]+
              T*(coeffShr[4]+
              T* coeffShr[5])));

          cp2=   coeffShr[8 ]+
              T*(coeffShr[9 ]+
              T*(coeffShr[10]+
              T*(coeffShr[11]+
              T* coeffShr[12])));

          if(T < Tmid)
          {
              cp2 = cp1;
          }
          Cp_R_dev[reactorid+speciesid*nReactors] = cp2;
      }
  }
}



void nasa_poly_group_getG_RT_CUDA(const double T, double * G_RT_dev, const int nGroupSpc, const double * thermoCoeff_dev, cudaStream_t stream)
{
    int nThreads = THREADS_PER_BLOCK;
    int nBlocks = (nGroupSpc+nThreads-1)/nThreads;
    cuda_get_G_RT<<<nBlocks, nThreads, 0, stream>>> (nGroupSpc, T, thermoCoeff_dev, G_RT_dev );
#ifdef ZERORK_FULL_DEBUG
    cudaDeviceSynchronize();
    checkCudaError(cudaGetLastError(),"cuda_get_G_RT");
#endif
}

void nasa_poly_group_getG_RT_CUDA_mr(const int nReactors, const double *T_dev, double *G_RT_dev, const int nGroupSpc, const double * thermoCoeff_dev, cudaStream_t stream)
{
    int nThreads = std::min(nReactors,MAX_THREADS_PER_BLOCK);
    dim3 nBlocks((nReactors+nThreads-1)/nThreads,nGroupSpc);
    cuda_get_G_RT_mr<<<nBlocks, nThreads, 0, stream>>>(nReactors, nGroupSpc, T_dev, thermoCoeff_dev, G_RT_dev);
#ifdef ZERORK_FULL_DEBUG
    cudaDeviceSynchronize();
    checkCudaError(cudaGetLastError(),"cuda_get_G_RT_mr");
#endif
}

void nasa_poly_group_getH_RT_CUDA_mr(const int nReactors, const double *T_dev, double *H_RT_dev, const int nGroupSpc, const double * thermoCoeff_dev, cudaStream_t stream)
{
    int nThreads = std::min(nReactors,MAX_THREADS_PER_BLOCK);
    dim3 nBlocks((nReactors+nThreads-1)/nThreads,nGroupSpc);
    cuda_get_H_RT_mr<<<nBlocks, nThreads, 0, stream>>>(nReactors, nGroupSpc, T_dev, thermoCoeff_dev, H_RT_dev);
#ifdef ZERORK_FULL_DEBUG
    cudaDeviceSynchronize();
    checkCudaError(cudaGetLastError(),"cuda_get_H_RT_mr");
#endif
}

void nasa_poly_group_getCp_R_CUDA_mr(const int nReactors, const double *T_dev, double *Cp_R_dev, const int nGroupSpc, const double * thermoCoeff_dev, cudaStream_t stream)
{
    int nThreads = std::min(nReactors,MAX_THREADS_PER_BLOCK);
    dim3 nBlocks((nReactors+nThreads-1)/nThreads,nGroupSpc);
    cuda_get_Cp_R_mr<<<nBlocks, nThreads, 0, stream>>>(nReactors, nGroupSpc, T_dev, thermoCoeff_dev, Cp_R_dev);
#ifdef ZERORK_FULL_DEBUG
    cudaDeviceSynchronize();
    checkCudaError(cudaGetLastError(),"cuda_get_Cp_R_mr");
#endif
}

} // namespace zerork
