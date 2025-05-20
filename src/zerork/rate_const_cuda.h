#ifndef ZERORK_RATE_CONST_CUDA_H
#define ZERORK_RATE_CONST_CUDA_H

#include "nasa_poly_cuda.h"
#include "rate_const.h"
#include <cuda_runtime.h>


namespace zerork {

class rate_const_cuda : public rate_const
{
 public:
  rate_const_cuda(ckr::CKReader *ckrobj, info_net *netobj, nasa_poly_group *tobj, int nReactorsMax = 64);
  virtual ~rate_const_cuda();

  void updateK_CUDA_mr(const int nReactors, const double T_dev[], const double * C_dev, double * K_dev);

  int nReactorsMax() { return m_nReactorsMax; };

 private:
  double *logAfact_dev;
  double *Tpow_dev;
  double *Tact_dev;
  double *Tmulti_dev;
  double *Csum_dev;

  //Equilibrium rates on GPU
  double *Gibbs_RT_dev;
  int keqPadSize;
  int *fromKeq_reacIdx_dev;
  int *fromKeq_prodIdx_dev;
  int *fromKeq_stepIdx_dev;
  int *fromKeq_fwdStepIdx_dev;
  double *fromKeq_nDelta_dev;
  double *fromKeq_stoich_fwd_dev;
  double *fromKeq_stoich_rev_dev;

  //Third body rates on GPU
  int maxThirdBodySpc; //max of both third body-only and falloff rxns
  int *thirdBody_nEnhanced_dev;
  int *thirdBody_etbSpcIdx_dev;
  int *thirdBody_fwdStepIdx_dev;
  int *thirdBody_revStepIdx_dev;
  double *thirdBody_etbSpcEff_dev;

  //Falloff rates on GPU
  int maxFalloffParams;
  int *falloff_falloffType_dev;
  int *falloff_nEnhanced_dev;
  int *falloff_falloffSpcIdx_dev;
  int *falloff_etbSpcIdx_dev;
  int *falloff_fwdStepIdx_dev;
  int *falloff_revStepIdx_dev;
  double *falloff_etbSpcEff_dev;
  double *falloff_param_dev;

  //PLog rates on GPU
  int nPLogInterpolationStepRev;
  int PLog_maxPressurePoints;
  double *PLog_signAfact_dev;
  double *PLog_logAfact_dev;
  double *PLog_Tpow_dev;
  double *PLog_Tact_dev;
  double *PLog_pressure_points_dev;
  int *PLog_stepIdxs_dev;

  //Streams
  cudaStream_t arrhStream,kEqStream,thirdBodyStream,falloffStream,PLogStream;

  void setGpuParams();

  int m_nReactorsMax;
};

} // namespace zerork

#endif
