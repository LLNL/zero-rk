
#ifndef ZERORK_RATE_CONST_KERNELS_H
#define ZERORK_RATE_CONST_KERNELS_H

namespace zerork {

void rate_const_updateArrheniusStep_CUDA_mr(const int nReactors,const int nStep,
    double *K_dev, const double *logAfact_dev, const double *Tpow_dev, const double *Tact_dev,
    const double *T_dev, cudaStream_t arrhStream);

void rate_const_concentrationSum(const int nReactors, const int nSpc, const double * C_dev, double *Csum_dev, cudaStream_t thirdBodyStream);

void rate_const_updateFromKeqStep_CUDA_mr(const int nReactors, const int nSpc, const int nStep,
          const int nFromKeqStep, const int keqPadSize, const int *fromKeq_reacIdx_dev,
          const int * fromKeq_prodIdx_dev, const int *fromKeq_stepIdx_dev,
          const double *fromKeq_nDelta_dev, const double *Gibbs_RT_dev, double *K_dev,
          const double *T_dev, cudaStream_t kEqStream);

void rate_const_updateFromKeqStepNI_CUDA_mr(const int nReactors, const int nSpc, const int nStep,
          const int nFromKeqStep, const int keqPadSize, const int *fromKeq_reacIdx_dev,
          const int * fromKeq_prodIdx_dev, const int *fromKeq_stepIdx_dev,
          const double *fromKeq_nDelta_dev, const double *fromKeq_stoich_fwd_dev,
          const double *fromKeq_stoich_rev_dev, const double *Gibbs_RT_dev, double *K_dev,
          const double *T_dev, cudaStream_t kEqStream);

void rate_const_updateThirdBody_CUDA_mr(const int nReactors, const int nSpc, const int nStep, 
    const int nThirdBodyRxn, const int maxThirdBodySpc, const int *thirdBody_nEnhanced_dev,
    const int *thirdBody_etbSpcIdx_dev, const int *thirdBody_fwdStepIdx_dev,
    const int *thirdBody_revStepIdx_dev, const double *thirdBody_etbSpcEff_dev, 
    const double *C_dev, const double *Csum_dev, double *K_dev, cudaStream_t thirdBodyStream);

void rate_const_updateFalloff_CUDA_mr(const int nReactors, const int nSpc, const int nStep,
         const int nFalloffRxn, const int maxThirdBodySpc, const int maxFalloffParams,
         const int *falloff_falloffType_dev, const int *falloff_nEnhanced_dev,
         const int *falloff_falloffSpcIdx_dev,
         const int *falloff_etbSpcIdx_dev, const int *falloff_fwdStepIdx_dev,
         const int *falloff_revStepIdx_dev, const double *falloff_etbSpcEff_dev,
         const double *falloff_param_dev, const double *logAfact_dev,
         const double *Tpow_dev, const double *Tact_dev, const double *C_dev,
         const double *Csum_dev, const double *T_dev, double *K_dev,
         cudaStream_t falloffStream);

void rate_const_updatePLogInterpolationStep_CUDA_mr(
  const int nReactors, const int nPLogInterpolationStep,
  const int PLog_maxPressureRanges,
  const int *PLog_stepIdxs_dev, const double *pressure_points_dev,
  const double * plog_signAfact_dev,
  const double * plog_logAfact_dev, const double *plog_Tpow_dev,
  const double * plog_Tact_dev, const double *Csum_dev,
  const double *T_dev, double *K_dev, cudaStream_t PLogStream);

} // namespace zerork

#endif
