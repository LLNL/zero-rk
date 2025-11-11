#include <algorithm> //std::min
#include "rate_const.h"
#include "rate_const_kernels.h"
#include "zerork_cuda_defs.h"
#include "constants.h"

namespace zerork {

void __global__ updateArrheniusStep_CUDA_mr
(
    const int nReactors,
    const int nStep,
    double *K_dev,
    const double *logAfact_dev,
    const double *Tpow_dev,
    const double *Tact_dev,
    const double *Tmulti_dev
)
{
    int reactorid = blockIdx.x*blockDim.x + threadIdx.x;
    int stepid = blockIdx.y*blockDim.y + threadIdx.y;

    extern __shared__ double shrMem[];
    double *logAfact_th= &shrMem[0];
    double *Tpow_th = &shrMem[blockDim.y];
    double *Tact_th = &shrMem[2*blockDim.y];
    double invT;
    double log_e_T;

    if(stepid < nStep)
    {
        if(reactorid < nReactors)
        {
            if(threadIdx.x == 0)
            {
                logAfact_th[threadIdx.y] = logAfact_dev[stepid];
                Tpow_th[threadIdx.y] = Tpow_dev[stepid];
                Tact_th[threadIdx.y] = Tact_dev[stepid];
            }
	}
    }
    __syncthreads();
    if(stepid < nStep)
    {
        if(reactorid < nReactors)
        {
            double T = Tmulti_dev[reactorid];
            invT= 1/T;
            log_e_T = log(T);

            K_dev[reactorid+nReactors*stepid] = exp( logAfact_th[threadIdx.y] 
                                  + Tpow_th[threadIdx.y]*log_e_T
                                  - Tact_th[threadIdx.y]*invT );
        }
    }
}


void __global__ updateFromKeqStep_CUDA_mr
(
    const int nReactors,
    const int nSpc,
    const int nStep,
    const int nFromKeqStep,
    const int nSpcPerStep,
    const int *fromKeq_reacIdx_dev,
    const int *fromKeq_prodIdx_dev,
    const int *fromKeq_stepIdx_dev,
    const double *fromKeq_nDelta_dev,
    const double *Gibbs_RT_dev,
    double *K_dev,
    const double *Tmulti_dev 
)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int reactorid = blockIdx.y*blockDim.y + threadIdx.y;
  extern int __shared__ reac_prodIdxs[];
  int *reacIdxs = &reac_prodIdxs[0];
  int *prodIdxs = &reac_prodIdxs[blockDim.x*nSpcPerStep];
  if(tid < nFromKeqStep)
  {
      if(reactorid < nReactors)
      {
          if(threadIdx.y == 0)
          {
              for(int j=0; j<nSpcPerStep; ++j )
              {
                  reacIdxs[threadIdx.x*nSpcPerStep + j] = fromKeq_reacIdx_dev[tid*nSpcPerStep+j];
                  prodIdxs[threadIdx.x*nSpcPerStep + j] = fromKeq_prodIdx_dev[tid*nSpcPerStep+j];
	      }
          }
      }
  }
  __syncthreads();
  if(tid < nFromKeqStep)
  {
      if(reactorid < nReactors)
      {
          double log_e_PatmInvRuT = log(P_ATM/(NIST_RU*Tmulti_dev[reactorid]));
          double gibbs_net = 0.0;

          for(int j=0; j<nSpcPerStep; ++j)
          {
              int currReacIdx = reacIdxs[threadIdx.x*nSpcPerStep+j]*nReactors+reactorid;
              gibbs_net -= Gibbs_RT_dev[currReacIdx];
              int currProdIdx = prodIdxs[threadIdx.x*nSpcPerStep+j]*nReactors+reactorid;
              gibbs_net += Gibbs_RT_dev[currProdIdx];
          }

          K_dev[reactorid + nReactors*fromKeq_stepIdx_dev[tid]] *= exp(gibbs_net-fromKeq_nDelta_dev[tid]*log_e_PatmInvRuT);
      }
  }
}

void __global__ updateFromKeqStepNI_CUDA_mr
(
    const int nReactors,
    const int nSpc,
    const int nStep,
    const int nFromKeqStep,
    const int nSpcPerStep,
    const int *fromKeq_reacIdx_dev,
    const int *fromKeq_prodIdx_dev,
    const int *fromKeq_stepIdx_dev,
    const double *fromKeq_nDelta_dev,
    const double *stoich_fwd_dev,
    const double *stoich_rev_dev,
    const double *Gibbs_RT_dev,
    double *K_dev,
    const double *Tmulti_dev 
)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int reactorid = blockIdx.y*blockDim.y + threadIdx.y;
  extern int __shared__ reac_prodIdxs[];
  int *reacIdxs = &reac_prodIdxs[0];
  int *prodIdxs = &reac_prodIdxs[blockDim.x*nSpcPerStep];
  if(tid < nFromKeqStep)
  {
      if(reactorid < nReactors)
      {
          //TODO: Consider packing stoich as well.  Need more shared mem for that
          if(threadIdx.y == 0)
          {
              for(int j=0; j<nSpcPerStep; ++j )
              {
                  reacIdxs[threadIdx.x*nSpcPerStep + j] = fromKeq_reacIdx_dev[tid*nSpcPerStep+j];
                  prodIdxs[threadIdx.x*nSpcPerStep + j] = fromKeq_prodIdx_dev[tid*nSpcPerStep+j];
	      }
          }
      }
  }
  __syncthreads();
  if(tid < nFromKeqStep)
  {
      if(reactorid < nReactors)
      {
          double log_e_PatmInvRuT = log(P_ATM/(NIST_RU*Tmulti_dev[reactorid]));
          double gibbs_net = 0.0;

          for(int j=0; j<nSpcPerStep; ++j)
          {
              const int arrIdx = threadIdx.x*nSpcPerStep+j;
              const int currReacIdx = reacIdxs[arrIdx]*nReactors+reactorid;
              gibbs_net -= stoich_fwd_dev[tid*nSpcPerStep+j]*Gibbs_RT_dev[currReacIdx];
              const int currProdIdx = prodIdxs[arrIdx]*nReactors+reactorid;
              gibbs_net += stoich_rev_dev[tid*nSpcPerStep+j]*Gibbs_RT_dev[currProdIdx];
          }

          K_dev[reactorid + nReactors*fromKeq_stepIdx_dev[tid]] *= exp(gibbs_net-fromKeq_nDelta_dev[tid]*log_e_PatmInvRuT);
      }
  }
}

void __global__ updateThirdBody_CUDA_mr
(
    const int nReactors,
    const int nSpc,
    const int nStep,
    const int nThirdBodyRxn,
    const int nSpcPerStep,
    const int *thirdBody_nEnhanced_dev,
    const int *thirdBody_etbSpcIdx_dev,
    const int *thirdBody_fwdStepIdx_dev,
    const int *thirdBody_revStepIdx_dev,
    const double *thirdBody_etbSpcEff_dev,
    const double *C_dev,
    const double *Csum_dev,
    double *K_dev
)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int reactorid = blockIdx.y*blockDim.y + threadIdx.y;
  if(tid < nThirdBodyRxn)
  {
      if(reactorid < nReactors)
      {
          double Cmult=Csum_dev[reactorid];
          for(int k=0; k<thirdBody_nEnhanced_dev[tid]; k++)
          {
              Cmult += C_dev[thirdBody_etbSpcIdx_dev[tid*nSpcPerStep + k]*nReactors+reactorid]*
                             thirdBody_etbSpcEff_dev[tid*nSpcPerStep + k];
	  }
          Cmult = Cmult;
          K_dev[thirdBody_fwdStepIdx_dev[tid]*nReactors+reactorid] *= Cmult;
          if(thirdBody_revStepIdx_dev[tid] >= 0)
          {
              K_dev[thirdBody_revStepIdx_dev[tid]*nReactors+reactorid] *= Cmult;
          }
      }
  }
}


void __global__ updateFalloff_CUDA_mr
(
    const int nReactors,
    const int nSpc,
    const int nStep,
    const int nFalloffRxn,
    const int nSpcPerStep,
    const int nFalloffParams,
    const int *falloff_falloffType_dev,
    const int *falloff_nEnhanced_dev,
    const int *falloff_falloffSpcIdx_dev,
    const int *falloff_etbSpcIdx_dev,
    const int *falloff_fwdStepIdx_dev,
    const int *falloff_revStepIdx_dev,
    const double *falloff_etbSpcEff_dev,
    const double *falloff_param_dev,
    const double *logAfact_dev,
    const double *Tpow_dev,
    const double *Tact_dev,
    const double *C_dev,
    const double *Csum_dev,
    const double *Tmulti_dev,
    double *K_dev
)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  int reactorid = blockIdx.y*blockDim.y + threadIdx.y;
  if(tid < nFalloffRxn)
  {
      if(reactorid < nReactors)
      {
          int fwdArrhIdx;
          double Pr,Klow,log_10_Pr,fTerm,Fcenter,Tcurrent,Cmult,log_e_Tcurrent,invTcurrent,
                 Kfwd,nTerm,Pcorr;
          Tcurrent = Tmulti_dev[reactorid];
          log_e_Tcurrent = log(Tcurrent);
          invTcurrent = 1.0/Tcurrent;

          if(falloff_falloffSpcIdx_dev[tid]>=0)
          {
             Cmult=C_dev[falloff_falloffSpcIdx_dev[tid]*nReactors+reactorid];
          }
          else
          {
              Cmult=Csum_dev[reactorid];
              for(int k=0; k<falloff_nEnhanced_dev[tid]; k++)
              {
                  Cmult+=C_dev[falloff_etbSpcIdx_dev[tid*nSpcPerStep+k]*nReactors+reactorid]*
                               falloff_etbSpcEff_dev[tid*nSpcPerStep+k];
              }
          }

          Klow = exp(falloff_param_dev[tid*nFalloffParams+0]
                    +falloff_param_dev[tid*nFalloffParams+1]*log_e_Tcurrent
                    -falloff_param_dev[tid*nFalloffParams+2]*invTcurrent);
//      Pr = Klow*Cmult/Kwork[falloffRxnList[j].fwdStepIdx];
          fwdArrhIdx = falloff_fwdStepIdx_dev[tid];
          Kfwd = exp(  logAfact_dev[fwdArrhIdx] 
                        + Tpow_dev[fwdArrhIdx]*log_e_Tcurrent
                        - Tact_dev[fwdArrhIdx]*invTcurrent );
          Pr = Klow*Cmult/Kfwd;
          if(Pr < 1.0e-300) {Pr = 1.0e-300;} // ck SMALL constant
          log_10_Pr=log10(Pr);

          fTerm=1.0; // default is Lindemann
          if(falloff_falloffType_dev[tid]==TROE_FOUR_PARAMS)
          {
              // Troe - 4 parameter fit
              //Fcenter=(1.0-falloff_param_dev[tid*nFalloffParams+3])
              //          *exp(-Tcurrent/falloff_param_dev[tid*nFalloffParams+4])
              //       +falloff_param_dev[tid*nFalloffParams+3]
              //          *exp(-Tcurrent/falloff_param_dev[tid*nFalloffParams+5])
              //       +exp(-falloff_param_dev[tid*nFalloffParams+6]*invTcurrent);
              Fcenter = 0.0;
              if(falloff_param_dev[tid*nFalloffParams+4]!=0) {
                Fcenter += (1.0-falloff_param_dev[tid*nFalloffParams+3])
                          *exp(-Tcurrent/falloff_param_dev[tid*nFalloffParams+4]);
              }
              if(falloff_param_dev[tid*nFalloffParams+5]!=0) {
                Fcenter += falloff_param_dev[tid*nFalloffParams+3]
                          *exp(-Tcurrent/falloff_param_dev[tid*nFalloffParams+5]);
              }
              Fcenter += exp(-falloff_param_dev[tid*nFalloffParams+6]*invTcurrent);

              if(Fcenter < 1.0e-300)
                {Fcenter=1.0e-300;}
              fTerm=log10(Fcenter);
              nTerm=0.75-1.27*fTerm;
              log_10_Pr-=(0.4+0.67*fTerm);
              log_10_Pr=log_10_Pr/(nTerm-0.14*log_10_Pr);
              log_10_Pr*=log_10_Pr;
              fTerm/=(1.0+log_10_Pr);
              fTerm=pow(10.,fTerm);
          } else if(falloff_falloffType_dev[tid]==TROE_THREE_PARAMS) {
              // Troe - 3 parameter fit
              //Fcenter=(1.0-falloff_param_dev[tid*nFalloffParams+3])
              //          *exp(-Tcurrent/falloff_param_dev[tid*nFalloffParams+4])
              //       +falloff_param_dev[tid*nFalloffParams+3]
              //          *exp(-Tcurrent/falloff_param_dev[tid*nFalloffParams+5]);
              Fcenter = 0.0;
              if(falloff_param_dev[tid*nFalloffParams+4]!=0) {
                Fcenter += (1.0-falloff_param_dev[tid*nFalloffParams+3])
                          *exp(-Tcurrent/falloff_param_dev[tid*nFalloffParams+4]);
              }
              if(falloff_param_dev[tid*nFalloffParams+5]!=0) {
                Fcenter += falloff_param_dev[tid*nFalloffParams+3]
                          *exp(-Tcurrent/falloff_param_dev[tid*nFalloffParams+5]);
              }

              if(Fcenter < 1.0e-300)
                {Fcenter=1.0e-300;}
              fTerm=log10(Fcenter);
              nTerm=0.75-1.27*fTerm;
              log_10_Pr-=(0.4+0.67*fTerm);
              log_10_Pr=log_10_Pr/(nTerm-0.14*log_10_Pr);
              log_10_Pr*=log_10_Pr;
              fTerm/=(1.0+log_10_Pr);
              fTerm=pow(10.,fTerm);
          }

          Pcorr = fTerm*Pr/(1.0+Pr);

          K_dev[falloff_fwdStepIdx_dev[tid]*nReactors+reactorid] *= Pcorr;
          if(falloff_revStepIdx_dev[tid] >= 0)
          {
              K_dev[falloff_revStepIdx_dev[tid]*nReactors+reactorid] *= Pcorr;
          }
      }
  }
}


void __global__  updatePLogInterpolationStep_CUDA_mr
(
  const int nReactors,
  const int nPLogInterpolationStep,
  const int maxpoints,
  const int *PLog_stepIdxs_dev, 
  const double *pressure_points,
  const double *sign_Afactor,
  const double *log_e_Afactor,
  const double *temperature_power,
  const double *activation_temperature,
  const double *Csum_dev,
  const double *T_dev,
  double *K_dev
)
{
    int stepid = blockIdx.x*blockDim.x + threadIdx.x;
    int reactorid = blockIdx.y*blockDim.y + threadIdx.y;
    if(stepid < nPLogInterpolationStep)
    {
        if(reactorid < nReactors)
        {
            double Treactor = T_dev[reactorid];
            double Preactor = NIST_RU*Csum_dev[reactorid]*Treactor;
            double log_e_pressure = log(Preactor);
            double log_e_temperature = log(Treactor);
            double inv_temperature = 1.0/Treactor;

            int rangeid = 0;
            //TODO : load pressure points into shrMem
            //N.B. : pressure points are set to reflect extrapolation choice in rate_const_cuda.cpp
            while(  rangeid+1 < maxpoints-1 //+1 so we don't increment again, -1 because range not points
                  && Preactor > pressure_points[stepid*maxpoints+rangeid+1]
                  &&   0.0    <= pressure_points[stepid*maxpoints+rangeid+2]) {++rangeid;}
            while( rangeid > 0 && pressure_points[stepid*maxpoints+rangeid] == pressure_points[stepid*maxpoints+rangeid-1] )
                {--rangeid;}

            double P_pt1 = pressure_points[stepid*maxpoints+rangeid];
            double K_pt1 = sign_Afactor[stepid*maxpoints+rangeid]*exp(log_e_Afactor[stepid*maxpoints+rangeid] + 
                                          temperature_power[stepid*maxpoints+rangeid]*log_e_temperature -
                                          activation_temperature[stepid*maxpoints+rangeid]*inv_temperature);

            while(  rangeid+1 < maxpoints-1
                  && pressure_points[stepid*maxpoints+rangeid] == pressure_points[stepid*maxpoints+rangeid+1])
            {
               ++rangeid;
               K_pt1 += sign_Afactor[stepid*maxpoints+rangeid]*exp(log_e_Afactor[stepid*maxpoints+rangeid] + 
                        temperature_power[stepid*maxpoints+rangeid]*log_e_temperature -
                        activation_temperature[stepid*maxpoints+rangeid]*inv_temperature);
            }

            double P_pt2 = pressure_points[stepid*maxpoints+rangeid+1];
            double K_pt2 = sign_Afactor[stepid*maxpoints+rangeid+1]*exp(log_e_Afactor[stepid*maxpoints+rangeid+1] + 
               temperature_power[stepid*maxpoints+rangeid+1]*log_e_temperature -
               activation_temperature[stepid*maxpoints+rangeid+1]*inv_temperature);

            while(  rangeid+2 < maxpoints-1 
                  && pressure_points[stepid*maxpoints+rangeid+1] == pressure_points[stepid*maxpoints+rangeid+2])
            {
               ++rangeid;
               K_pt2 += sign_Afactor[stepid*maxpoints+rangeid+1]*exp(log_e_Afactor[stepid*maxpoints+rangeid+1] + 
                        temperature_power[stepid*maxpoints+rangeid+1]*log_e_temperature -
                        activation_temperature[stepid*maxpoints+rangeid+1]*inv_temperature);
            }

            double interpolation_exponent = log(K_pt2/K_pt1) / log(P_pt2/P_pt1);
            K_dev[PLog_stepIdxs_dev[stepid]*nReactors+reactorid] = K_pt1*pow(Preactor/P_pt1,interpolation_exponent);
      }
  }
}


void __global__ concentrationSum
(
    int nReactors,
    int nSpc,
    const double *C_dev,
    double *Csum_dev
)
{
  int reactorid = blockIdx.x*blockDim.x + threadIdx.x;
  if(reactorid < nReactors)
  {
      double Csum_local = 0.0;
      for(int j = 0; j < nSpc; ++j)
      {
          Csum_local += C_dev[nReactors*j+reactorid];
      }
      Csum_dev[reactorid] = Csum_local;
  }
}


void rate_const_updateArrheniusStep_CUDA_mr(const int nReactors,const int nStep,
    double *K_dev, const double *logAfact_dev, const double *Tpow_dev, const double *Tact_dev,
    const double *T_dev, cudaStream_t arrhStream)
{
  //Tuned on P100
  int threadsX = 256;
  int threadsY = 1;
  dim3 nThreads2D(threadsX,threadsY);

  int nBlocksX = (nReactors+threadsX-1)/threadsX;
  int nBlocksY = (nStep+threadsY-1)/threadsY;
  dim3 nBlocks2D(nBlocksX,nBlocksY);

  size_t shrMemSize = (3*nThreads2D.y)*sizeof(double);
  updateArrheniusStep_CUDA_mr<<<nBlocks2D,nThreads2D, shrMemSize,arrhStream>>>
                 (nReactors,nStep,K_dev,logAfact_dev,Tpow_dev,Tact_dev,T_dev);
#ifdef ZERORK_FULL_DEBUG
  cudaDeviceSynchronize();
  checkCudaError(cudaGetLastError(),"updateArrheniusStepCUDA_mr");
#endif
}

void rate_const_concentrationSum(const int nReactors, const int nSpc, const double * C_dev, double *Csum_dev, cudaStream_t thirdBodyStream)
{
  int nThreads = std::min(THREADS_PER_BLOCK,nReactors);
  int nBlocks = (nReactors+nThreads-1)/nThreads;
  concentrationSum<<<nBlocks,nThreads,0,thirdBodyStream>>>(nReactors, nSpc, C_dev, Csum_dev);
#ifdef ZERORK_FULL_DEBUG
  cudaDeviceSynchronize();
  checkCudaError(cudaGetLastError(),"concentrationSum");
#endif
}


void rate_const_updateFromKeqStep_CUDA_mr(const int nReactors, const int nSpc, const int nStep,
          const int nFromKeqStep, const int keqPadSize, const int *fromKeq_reacIdx_dev,
          const int * fromKeq_prodIdx_dev, const int *fromKeq_stepIdx_dev,
          const double *fromKeq_nDelta_dev, const double *Gibbs_RT_dev, double *K_dev,
          const double *T_dev, cudaStream_t kEqStream)
{
  if( nFromKeqStep == 0 ) return;
  int threadsX = 1;
  int threadsY = std::min(nReactors,MAX_THREADS_PER_BLOCK/threadsX);
  dim3 nThreads2D(threadsX,threadsY);
  dim3 nBlocks2D;
  nBlocks2D.x = (nFromKeqStep+threadsX-1)/threadsX;
  nBlocks2D.y = (nReactors+threadsY-1)/threadsY;
  updateFromKeqStep_CUDA_mr<<<nBlocks2D,nThreads2D, nThreads2D.x*keqPadSize*2*sizeof(int),kEqStream>>>
  (
      nReactors,
      nSpc,
      nStep,
      nFromKeqStep,
      keqPadSize,
      fromKeq_reacIdx_dev,
      fromKeq_prodIdx_dev,
      fromKeq_stepIdx_dev,
      fromKeq_nDelta_dev,
      Gibbs_RT_dev,
      K_dev,
      T_dev
  );
#ifdef ZERORK_FULL_DEBUG
  cudaDeviceSynchronize();
  checkCudaError(cudaGetLastError(),"updateFromKeq_CUDA_mr");
#endif
}

void rate_const_updateFromKeqStepNI_CUDA_mr(const int nReactors, const int nSpc, const int nStep,
          const int nFromKeqStep, const int keqPadSize, const int *fromKeq_reacIdx_dev,
          const int * fromKeq_prodIdx_dev, const int *fromKeq_stepIdx_dev,
          const double *fromKeq_nDelta_dev, const double *fromKeq_stoich_fwd_dev,
          const double *fromKeq_stoich_rev_dev, const double *Gibbs_RT_dev, double *K_dev,
          const double *T_dev, cudaStream_t kEqStream)
{
  if( nFromKeqStep == 0 ) return;
  int threadsX = 1;
  int threadsY = std::min(nReactors,MAX_THREADS_PER_BLOCK/threadsX);
  dim3 nThreads2D(threadsX,threadsY);
  dim3 nBlocks2D;
  nBlocks2D.x = (nFromKeqStep+threadsX-1)/threadsX;
  nBlocks2D.y = (nReactors+threadsY-1)/threadsY;
  updateFromKeqStepNI_CUDA_mr<<<nBlocks2D,nThreads2D, nThreads2D.x*keqPadSize*2*sizeof(int),kEqStream>>>
  (
      nReactors,
      nSpc,
      nStep,
      nFromKeqStep,
      keqPadSize,
      fromKeq_reacIdx_dev,
      fromKeq_prodIdx_dev,
      fromKeq_stepIdx_dev,
      fromKeq_nDelta_dev,
      fromKeq_stoich_fwd_dev,
      fromKeq_stoich_rev_dev,
      Gibbs_RT_dev,
      K_dev,
      T_dev
  );
#ifdef ZERORK_FULL_DEBUG
  cudaDeviceSynchronize();
  checkCudaError(cudaGetLastError(),"updateFromKeq_CUDA_mr");
#endif
}


void rate_const_updateThirdBody_CUDA_mr(const int nReactors, const int nSpc, const int nStep, 
    const int nThirdBodyRxn, const int maxThirdBodySpc, const int *thirdBody_nEnhanced_dev,
    const int *thirdBody_etbSpcIdx_dev, const int *thirdBody_fwdStepIdx_dev,
    const int *thirdBody_revStepIdx_dev, const double *thirdBody_etbSpcEff_dev, 
    const double *C_dev, const double *Csum_dev, double *K_dev, cudaStream_t thirdBodyStream)
{
  if( nThirdBodyRxn == 0 ) return;
  int threadsX = 1;
  int threadsY = std::min(nReactors,MAX_THREADS_PER_BLOCK/threadsX);
  dim3 nThreads2D(threadsX,threadsY);
  dim3 nBlocks2D;
  nBlocks2D.x = (nThirdBodyRxn+threadsX-1)/threadsX;
  nBlocks2D.y = (nReactors+threadsY-1)/threadsY;
  updateThirdBody_CUDA_mr<<<nBlocks2D,nThreads2D,0,thirdBodyStream>>>
  (
     nReactors,
     nSpc,
     nStep,
     nThirdBodyRxn,
     maxThirdBodySpc,
     thirdBody_nEnhanced_dev,
     thirdBody_etbSpcIdx_dev,
     thirdBody_fwdStepIdx_dev,
     thirdBody_revStepIdx_dev,
     thirdBody_etbSpcEff_dev,
     C_dev,
     Csum_dev,
     K_dev
  );
#ifdef ZERORK_FULL_DEBUG
  cudaDeviceSynchronize();
  checkCudaError(cudaGetLastError(),"updateThirdBody_CUDA_mr");
#endif
}

void rate_const_updateFalloff_CUDA_mr(const int nReactors, const int nSpc, const int nStep,
         const int nFalloffRxn, const int maxThirdBodySpc, const int maxFalloffParams,
         const int *falloff_falloffType_dev, const int *falloff_nEnhanced_dev,
         const int *falloff_falloffSpcIdx_dev,
         const int *falloff_etbSpcIdx_dev, const int *falloff_fwdStepIdx_dev,
         const int *falloff_revStepIdx_dev, const double *falloff_etbSpcEff_dev,
         const double *falloff_param_dev, const double *logAfact_dev,
         const double *Tpow_dev, const double *Tact_dev, const double *C_dev,
         const double *Csum_dev, const double *T_dev, double *K_dev,
         cudaStream_t falloffStream)
{
  if( nFalloffRxn == 0 ) return;
  int threadsX = 1;
  int threadsY = std::min(nReactors,MAX_THREADS_PER_BLOCK/threadsX);
  dim3 nThreads2D(threadsX,threadsY);
  dim3 nBlocks2D;
  nBlocks2D.x = (nFalloffRxn+threadsX-1)/threadsX;
  nBlocks2D.y = (nReactors+threadsY-1)/threadsY;
  updateFalloff_CUDA_mr<<<nBlocks2D,nThreads2D,0,falloffStream>>>
  (
     nReactors,
     nSpc,
     nStep,
     nFalloffRxn,
     maxThirdBodySpc,
     maxFalloffParams,
     falloff_falloffType_dev,
     falloff_nEnhanced_dev,
     falloff_falloffSpcIdx_dev,
     falloff_etbSpcIdx_dev,
     falloff_fwdStepIdx_dev,
     falloff_revStepIdx_dev,
     falloff_etbSpcEff_dev,
     falloff_param_dev,
     logAfact_dev,
     Tpow_dev,
     Tact_dev,
     C_dev,
     Csum_dev,
     T_dev,
     K_dev
  );
#ifdef ZERORK_FULL_DEBUG
  cudaDeviceSynchronize();
  checkCudaError(cudaGetLastError(),"updateFalloff_CUDA_mr");
#endif
}

void rate_const_updatePLogInterpolationStep_CUDA_mr(
  const int nReactors, const int nPLogInterpolationStep, 
  const int PLog_maxPressurePoints,
  const int *PLog_stepIdxs_dev, const double *pressure_points_dev,
  const double * plog_signAfact_dev,
  const double * plog_logAfact_dev, const double *plog_Tpow_dev,
  const double * plog_Tact_dev, const double *Csum_dev,
  const double *T_dev, double *K_dev, cudaStream_t PLogStream)
{
  if(nPLogInterpolationStep == 0) return;
  int threadsX = 1;
  int threadsY = std::min(nReactors,MAX_THREADS_PER_BLOCK/threadsX);
  dim3 nThreads2D(threadsX,threadsY);

  int nBlocksX = (nPLogInterpolationStep+threadsX-1)/threadsX;
  int nBlocksY = (nReactors+threadsY-1)/threadsY;
  dim3 nBlocks2D(nBlocksX,nBlocksY);

//  size_t shrMemSize = (3*nThreads2D.x + 2*nThreads2D.y)*sizeof(double);
//  size_t shrMemSize = (3*nThreads2D.x)*sizeof(double);
  size_t shrMemSize = 0;
  updatePLogInterpolationStep_CUDA_mr<<<nBlocks2D,nThreads2D, shrMemSize,PLogStream>>>
         (nReactors,  nPLogInterpolationStep, PLog_maxPressurePoints, PLog_stepIdxs_dev, 
          pressure_points_dev, plog_signAfact_dev, plog_logAfact_dev, plog_Tpow_dev, plog_Tact_dev,
          Csum_dev, T_dev, K_dev);
#ifdef ZERORK_FULL_DEBUG
  cudaDeviceSynchronize();
  checkCudaError(cudaGetLastError(),"updateArrheniusStepCUDA_mr");
#endif
}

} // namespace zerork
