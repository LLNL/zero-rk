#include <stdlib.h>
#include <math.h>
#include <limits>
#include "rate_const.h"
#include "rate_const_cuda.h"
#include "constants.h"
#include "fast_exps.h"
#include "zerork_cuda_defs.h"
#include "rate_const_kernels.h"

#define unlikely(expr) __builtin_expect(!!(expr), 0)
#define likely(expr) __builtin_expect(!!(expr), 1)

namespace zerork {

rate_const_cuda::rate_const_cuda(ckr::CKReader *ckrobj, info_net *netobj, 
			       nasa_poly_group *tobj, int nReactorsMax)
    : rate_const(ckrobj,netobj,tobj), m_nReactorsMax(nReactorsMax)
{
  delete [] Gibbs_RT; //Re-allocate to larger size
  Gibbs_RT = new double[nSpc*m_nReactorsMax]; //TODO: May be leaking this.  Not sure if rate_const() destructor is called.

  setGpuParams();
}

rate_const_cuda::~rate_const_cuda()
{
  hipFree(logAfact_dev);
  hipFree(Tpow_dev);
  hipFree(Tact_dev);
  hipFree(Tmulti_dev);
  hipFree(Csum_dev);
  hipFree(Gibbs_RT_dev);

  hipFree(fromKeq_reacIdx_dev);
  hipFree(fromKeq_prodIdx_dev);
  hipFree(fromKeq_stepIdx_dev);
  hipFree(fromKeq_nDelta_dev);
  hipFree(fromKeq_stoich_fwd_dev);
  hipFree(fromKeq_stoich_rev_dev);

  hipFree(thirdBody_fwdStepIdx_dev);
  hipFree(thirdBody_revStepIdx_dev);
  hipFree(thirdBody_nEnhanced_dev);
  hipFree(thirdBody_etbSpcIdx_dev);
  hipFree(thirdBody_etbSpcEff_dev);

  hipFree(falloff_falloffType_dev);
  hipFree(falloff_fwdStepIdx_dev);
  hipFree(falloff_revStepIdx_dev);
  hipFree(falloff_nEnhanced_dev);
  hipFree(falloff_falloffSpcIdx_dev);
  hipFree(falloff_etbSpcIdx_dev);
  hipFree(falloff_etbSpcEff_dev);
  hipFree(falloff_param_dev);

  hipFree(PLog_signAfact_dev);
  hipFree(PLog_logAfact_dev);
  hipFree(PLog_Tpow_dev);
  hipFree(PLog_Tact_dev);
  hipFree(PLog_pressure_points_dev);
  hipFree(PLog_stepIdxs_dev);

  hipStreamDestroy(arrhStream);
  hipStreamDestroy(kEqStream);
  hipStreamDestroy(thirdBodyStream);
  hipStreamDestroy(falloffStream);
  hipStreamDestroy(PLogStream);
}

void rate_const_cuda::setGpuParams()
{
  int j,k,m;

  //TODO: Convert C99 style array allocation to vectors
  double logAfact_tmp[nStep];
  double Tpow_tmp[nStep];
  double Tact_tmp[nStep];
  memset(logAfact_tmp,0,nStep*sizeof(double));
  memset(Tpow_tmp,0,nStep*sizeof(double));
  memset(Tact_tmp,0,nStep*sizeof(double));

  for(j=0; j<nArrheniusStep; ++j)
    {
      logAfact_tmp[arrheniusStepList[j].stepIdx] =
                   distinctArrheniusLogAfact[arrheniusStepList[j].arrheniusIdx];
      Tpow_tmp[arrheniusStepList[j].stepIdx] =
                   distinctArrheniusTpow[arrheniusStepList[j].arrheniusIdx];
      Tact_tmp[arrheniusStepList[j].stepIdx] =
                   distinctArrheniusTact[arrheniusStepList[j].arrheniusIdx];
    }

  keqPadSize = max(MAX_NREACS,MAX_NPRODS);
  // CPU memory to store fromKeqStepList values.
  int fromKeq_reacIdx_tmp[nFromKeqStep*keqPadSize];
  int fromKeq_prodIdx_tmp[nFromKeqStep*keqPadSize];
  int fromKeq_stepIdx_tmp[nFromKeqStep];
  double fromKeq_nDelta_tmp[nFromKeqStep];
  double fromKeq_stoich_fwd_tmp[nFromKeqStep*keqPadSize];
  double fromKeq_stoich_rev_tmp[nFromKeqStep*keqPadSize];
  for(j=0; j < nFromKeqStep; ++j)
  {
      const int fwdStepIdx = fromKeqStepList[j].fwdStepIdx;
      logAfact_tmp[fromKeqStepList[j].stepIdx]=logAfact_tmp[fwdStepIdx];
      Tpow_tmp[fromKeqStepList[j].stepIdx]=Tpow_tmp[fwdStepIdx];
      Tact_tmp[fromKeqStepList[j].stepIdx]=Tact_tmp[fwdStepIdx];

      fromKeq_stepIdx_tmp[j] = fromKeqStepList[j].stepIdx;
      fromKeq_nDelta_tmp[j] = fromKeqStepList[j].nDelta;
      if(!non_integer_network_.HasStep(fwdStepIdx)){ 
        for(k = 0; k < fromKeqStepList[j].nReac; ++k)
        {
            fromKeq_reacIdx_tmp[j*keqPadSize + k] = fromKeqStepList[j].reacSpcIdx[k];
            fromKeq_stoich_fwd_tmp[j*keqPadSize + k] =  1.0;
        }
        for(k = fromKeqStepList[j].nReac; k < keqPadSize; ++k)
        {
            fromKeq_reacIdx_tmp[j*keqPadSize + k] = nSpc;
            fromKeq_stoich_fwd_tmp[j*keqPadSize + k] =  0.0;
        }
        for(k = 0; k < fromKeqStepList[j].nProd; ++k)
        {
            fromKeq_prodIdx_tmp[j*keqPadSize + k] = fromKeqStepList[j].prodSpcIdx[k];
            fromKeq_stoich_rev_tmp[j*keqPadSize + k] =  1.0;
        }
        for(k = fromKeqStepList[j].nProd; k < keqPadSize; ++k)
        {
            fromKeq_prodIdx_tmp[j*keqPadSize + k] = nSpc;
            fromKeq_stoich_rev_tmp[j*keqPadSize + k] =  0.0;
        }
      } else {
        const int num_prod = non_integer_network_.GetNumProductsOfStep(fwdStepIdx);
        for(k = 0; k < num_prod; ++k) {
          fromKeq_prodIdx_tmp[j*keqPadSize + k] = non_integer_network_.GetProductIndexOfStep(fwdStepIdx,k);
          fromKeq_stoich_rev_tmp[j*keqPadSize + k] =  non_integer_network_.GetProductPowerOfStep(fwdStepIdx,k);
        }
        for(k = num_prod; k < keqPadSize; ++k) {
          fromKeq_prodIdx_tmp[j*keqPadSize + k] = nSpc;
          fromKeq_stoich_rev_tmp[j*keqPadSize + k] =  0.0;
        }
        const int num_reac = non_integer_network_.GetNumReactantsOfStep(fwdStepIdx);
        for(k = 0; k < num_reac; ++k) {
          fromKeq_reacIdx_tmp[j*keqPadSize + k] = non_integer_network_.GetReactantIndexOfStep(fwdStepIdx,k);
          fromKeq_stoich_fwd_tmp[j*keqPadSize + k] =  non_integer_network_.GetReactantPowerOfStep(fwdStepIdx,k);
        }
        for(k = num_reac; k < keqPadSize; ++k) {
          fromKeq_reacIdx_tmp[j*keqPadSize + k] = nSpc;
          fromKeq_stoich_fwd_tmp[j*keqPadSize + k] =  0.0;
        }
      }
  }

  //Third body rates on GPU
  maxThirdBodySpc = 0;
  for(j=0; j<nThirdBodyRxn; j++)
  {
      maxThirdBodySpc = max(maxThirdBodySpc, thirdBodyRxnList[j].nEnhanced);
  }
  for(j=0; j<nFalloffRxn; j++)
  {
      maxThirdBodySpc = max(maxThirdBodySpc, falloffRxnList[j].nEnhanced);
  }

  int thirdBody_fwdStepIdx_tmp[nThirdBodyRxn];
  int thirdBody_revStepIdx_tmp[nThirdBodyRxn];
  int thirdBody_nEnhanced_tmp[nThirdBodyRxn];
  int thirdBody_etbSpcIdx_tmp[nThirdBodyRxn*maxThirdBodySpc];
  double thirdBody_etbSpcEff_tmp[nThirdBodyRxn*maxThirdBodySpc];

  for(j=0; j<nThirdBodyRxn; j++)
  {
      thirdBody_fwdStepIdx_tmp[j] = thirdBodyRxnList[j].fwdStepIdx;
      thirdBody_revStepIdx_tmp[j] = thirdBodyRxnList[j].revStepIdx;
      thirdBody_nEnhanced_tmp[j] = thirdBodyRxnList[j].nEnhanced;
      for(k=0; k < thirdBody_nEnhanced_tmp[j]; ++k)
      {
          thirdBody_etbSpcIdx_tmp[j*maxThirdBodySpc+k] = thirdBodyRxnList[j].etbSpcIdx[k];
          thirdBody_etbSpcEff_tmp[j*maxThirdBodySpc+k] = thirdBodyRxnList[j].etbSpcEff[k];
      }
  }

  //Falloff rates on GPU
  maxFalloffParams = 7;
  int falloff_falloffType_tmp[nFalloffRxn];
  int falloff_fwdStepIdx_tmp[nFalloffRxn];
  int falloff_revStepIdx_tmp[nFalloffRxn];
  int falloff_nEnhanced_tmp[nFalloffRxn];
  int falloff_falloffSpcIdx_tmp[nFalloffRxn];
  int falloff_etbSpcIdx_tmp[nFalloffRxn*maxThirdBodySpc];
  double falloff_etbSpcEff_tmp[nFalloffRxn*maxThirdBodySpc];
  double falloff_param_tmp[nFalloffRxn*maxFalloffParams];
  memset(falloff_etbSpcIdx_tmp,0,nFalloffRxn*maxThirdBodySpc*sizeof(int));
  memset(falloff_etbSpcEff_tmp,0,nFalloffRxn*maxThirdBodySpc*sizeof(double));
  memset(falloff_param_tmp,0,nFalloffRxn*maxFalloffParams*sizeof(int));
  for(j=0; j<nFalloffRxn; j++)
  {
      falloff_falloffType_tmp[j] = falloffRxnList[j].falloffType;
      falloff_fwdStepIdx_tmp[j] = falloffRxnList[j].fwdStepIdx;
      falloff_revStepIdx_tmp[j] = falloffRxnList[j].revStepIdx;
      falloff_falloffSpcIdx_tmp[j] = falloffRxnList[j].falloffSpcIdx;
      falloff_nEnhanced_tmp[j] = falloffRxnList[j].nEnhanced;
      for(k=0; k < falloff_nEnhanced_tmp[j]; ++k)
      {
          falloff_etbSpcIdx_tmp[j*maxThirdBodySpc+k] = falloffRxnList[j].etbSpcIdx[k];
          falloff_etbSpcEff_tmp[j*maxThirdBodySpc+k] = falloffRxnList[j].etbSpcEff[k];
      }
      int nCurrParams = 3;
      if(falloff_falloffType_tmp[j] == TROE_THREE_PARAMS)
      {
          nCurrParams = 6;
      }
      else if(falloff_falloffType_tmp[j] == TROE_FOUR_PARAMS)
      {
          nCurrParams = 7;
      }
      for(k=0; k < nCurrParams; ++k)
      {
          falloff_param_tmp[j*maxFalloffParams+k] = falloffRxnList[j].param[k];
      }
  }

  //PLog reactions on GPU
  PLog_maxPressurePoints = 0;
  std::map<int,std::pair<int,int> > PLog_eq_step_map;
  for(j=0; j<nPLogInterpolationStep; ++j)
  {
    PLogReaction& rxn = plogInterpolationStepList[j];
    bool extrap = rxn.use_extrapolation();
    int p_points = 0;
    for(m=0; m<rxn.num_pressure_points(); ++m) {
       p_points += rxn.num_arrhenius_lines(m);
    }
    if(!extrap)
    {
       p_points += rxn.num_arrhenius_lines(0);
       p_points += rxn.num_arrhenius_lines(rxn.num_pressure_points()-1);
    }
    PLog_maxPressurePoints = std::max(p_points,PLog_maxPressurePoints);
    PLog_eq_step_map[rxn.step_index()] = std::make_pair(-1,j); //second pair are rev eq and PLog counter
  }

  nPLogInterpolationStepRev = 0;
  for(j=0; j<nFromKeqStep; ++j)
  {
    if(PLog_eq_step_map.count(fromKeqStepList[j].fwdStepIdx) == 1)
    {
      ++nPLogInterpolationStepRev;
      PLog_eq_step_map[fromKeqStepList[j].fwdStepIdx].first = fromKeqStepList[j].stepIdx;
    }
  }
  for(j=0; j<nPLogInterpolationStep; ++j)
  {
    PLogReaction& rxn = plogInterpolationStepList[j];
    if(PLog_eq_step_map[rxn.step_index()].first == -1) 
    {
      PLog_eq_step_map.erase(rxn.step_index());
    }
  }

  std::vector<double> plog_sign_Afact(PLog_maxPressurePoints*(nPLogInterpolationStep+nPLogInterpolationStepRev),0.0);
  std::vector<double> plog_logAfact(PLog_maxPressurePoints*(nPLogInterpolationStep+nPLogInterpolationStepRev),0.0);
  std::vector<double> plog_Tpow(PLog_maxPressurePoints*(nPLogInterpolationStep+nPLogInterpolationStepRev),0.0);
  std::vector<double> plog_Tact(PLog_maxPressurePoints*(nPLogInterpolationStep+nPLogInterpolationStepRev),0.0);
  //NOTE: Should be safe to assume no negative pressure points.  We use this is a 
  //      flag in the kernel to stop looping without having a separate piece of
  //      data we'd have to pull from global memory.  If we want to allow negative
  //      pressure points and/or negative pressures during simulation we will
  //      have to revisit this choice.
  std::vector<double> plog_pressure_points(PLog_maxPressurePoints*(nPLogInterpolationStep+nPLogInterpolationStepRev),-1.0);
  std::vector<int> plog_stepIdxs(nPLogInterpolationStep+nPLogInterpolationStepRev,-1);
  for(j=0; j<nPLogInterpolationStep; ++j)
  {
    PLogReaction& rxn = plogInterpolationStepList[j];
    plog_stepIdxs[j] = rxn.step_index();
    bool extrap = rxn.use_extrapolation();
    int idx=0;
    int line_idx = 0;
    if(!extrap)
    {
       for(m=0; m < rxn.num_arrhenius_lines(0); ++m) {
           plog_sign_Afact[PLog_maxPressurePoints*j+idx] = rxn.sign_afactor_[m];
           plog_logAfact[PLog_maxPressurePoints*j+idx] = rxn.log_e_afactor_[m];
           plog_Tpow[PLog_maxPressurePoints*j+idx] = rxn.temperature_power_[m];
           plog_Tact[PLog_maxPressurePoints*j+idx] = rxn.activation_temperature_[m];
           plog_pressure_points[PLog_maxPressurePoints*j+idx] = 0.0;
           ++idx;
           ++line_idx;
       }
    }
    line_idx=0;
    for(k=0; k < rxn.num_pressure_points(); ++k)
    {
       for(m=0; m < rxn.num_arrhenius_lines(k); ++m) {
         plog_sign_Afact[PLog_maxPressurePoints*j+idx] = rxn.sign_afactor_[line_idx];
         plog_logAfact[PLog_maxPressurePoints*j+idx] = rxn.log_e_afactor_[line_idx];
         plog_Tpow[PLog_maxPressurePoints*j+idx] = rxn.temperature_power_[line_idx];
         plog_Tact[PLog_maxPressurePoints*j+idx] = rxn.activation_temperature_[line_idx];
         plog_pressure_points[PLog_maxPressurePoints*j+idx] = rxn.pressure_points_[k];
         ++idx;
         ++line_idx;
       }
    }
    if(!extrap)
    {
       line_idx -= rxn.num_arrhenius_lines(k-1);
       for(m=0; m < rxn.num_arrhenius_lines(k-1); ++m) {
           plog_sign_Afact[PLog_maxPressurePoints*j+idx] = rxn.sign_afactor_[line_idx+m];
           plog_logAfact[PLog_maxPressurePoints*j+idx] = rxn.log_e_afactor_[line_idx+m];
           plog_Tpow[PLog_maxPressurePoints*j+idx] = rxn.temperature_power_[line_idx+m];
           plog_Tact[PLog_maxPressurePoints*j+idx] = rxn.activation_temperature_[line_idx+m];
           plog_pressure_points[PLog_maxPressurePoints*j+idx] = std::numeric_limits<double>::max();
           ++idx;
        }
    }
  }


  j=nPLogInterpolationStep;
  for(std::map<int,std::pair<int,int> >::iterator it = PLog_eq_step_map.begin();
                                  it != PLog_eq_step_map.end();
                                  ++it)
  {
    int revIdx = it->second.first;
    int jPLog = it->second.second;
    PLogReaction& rxn = plogInterpolationStepList[jPLog];
    plog_stepIdxs[j] = revIdx;

    bool extrap = rxn.use_extrapolation();
    int idx=0;
    int line_idx=0;
    if(!extrap)
    {
       for(m=0; m < rxn.num_arrhenius_lines(0); ++m) {
           plog_sign_Afact[PLog_maxPressurePoints*j+idx] = rxn.sign_afactor_[m];
           plog_logAfact[PLog_maxPressurePoints*j+idx] = rxn.log_e_afactor_[m];
           plog_Tpow[PLog_maxPressurePoints*j+idx] = rxn.temperature_power_[m];
           plog_Tact[PLog_maxPressurePoints*j+idx] = rxn.activation_temperature_[m];
           plog_pressure_points[PLog_maxPressurePoints*j+idx] = 0.0;
           ++idx;
           ++line_idx;
       }
    }
    line_idx = 0;
    for(k=0; k < rxn.num_pressure_points(); ++k)
    {
       for(m=0; m < rxn.num_arrhenius_lines(k); ++m) {
           plog_sign_Afact[PLog_maxPressurePoints*j+idx] = rxn.sign_afactor_[line_idx];
           plog_logAfact[PLog_maxPressurePoints*j+idx] = rxn.log_e_afactor_[line_idx];
           plog_Tpow[PLog_maxPressurePoints*j+idx] = rxn.temperature_power_[line_idx];
           plog_Tact[PLog_maxPressurePoints*j+idx] = rxn.activation_temperature_[line_idx];
           plog_pressure_points[PLog_maxPressurePoints*j+idx] = rxn.pressure_points_[k];
           ++idx;
           ++line_idx;
       }
    }
    if(!extrap)
    {
       line_idx -= rxn.num_arrhenius_lines(k-1);
       for(m=0; m < rxn.num_arrhenius_lines(k-1); ++m) {
           plog_sign_Afact[PLog_maxPressurePoints*j+idx] = rxn.sign_afactor_[line_idx+m];
           plog_logAfact[PLog_maxPressurePoints*j+idx] = rxn.log_e_afactor_[line_idx+m];
           plog_Tpow[PLog_maxPressurePoints*j+idx] = rxn.temperature_power_[line_idx+m];
           plog_Tact[PLog_maxPressurePoints*j+idx] = rxn.activation_temperature_[line_idx+m];
           plog_pressure_points[PLog_maxPressurePoints*j+idx] = std::numeric_limits<double>::max();
           ++idx;
        }
    }
    ++j;
  }

  checkCudaError
  (
     hipMalloc((void**)&Tmulti_dev,sizeof(double)*m_nReactorsMax),
     "hipMalloc(... Tmulti_dev ...)"
  );

  checkCudaError
  (
     hipMalloc((void**)&Csum_dev,sizeof(double)*m_nReactorsMax),
     "hipMalloc(... Csum_dev ...)"
  );

  // NB : nSpc+1 is for padded space so we can unroll kEq calculation.  CF. fromKeq_{reac,prod}Idx_{tmp,dev}
  checkCudaError
  (
     hipMalloc((void**)&Gibbs_RT_dev,sizeof(double)*(nSpc+1)*m_nReactorsMax),
     "hipMalloc(... Gibbs_RT_dev ...)"
  );
  hipMemset(&Gibbs_RT_dev[nSpc*m_nReactorsMax],0,sizeof(double)*m_nReactorsMax);

  checkCudaError
  (
     hipMalloc((void**)&logAfact_dev,sizeof(double)*nStep),
     "hipMalloc(... logAfact_dev ...)"
  );
  checkCudaError
  (
     hipMalloc((void**)&Tpow_dev,sizeof(double)*nStep),
     "hipMalloc(... Tpow_dev ...)"
  );
  checkCudaError
  (
     hipMalloc((void**)&Tact_dev,sizeof(double)*nStep),
     "hipMalloc(... Tact_dev ...)"
  );
  hipMemcpy(logAfact_dev,logAfact_tmp,nStep*sizeof(double),hipMemcpyHostToDevice);
  hipMemcpy(Tpow_dev,Tpow_tmp,nStep*sizeof(double),hipMemcpyHostToDevice);
  hipMemcpy(Tact_dev,Tact_tmp,nStep*sizeof(double),hipMemcpyHostToDevice);


  //Keq on GPU
  checkCudaError
  (
     hipMalloc((void**)&fromKeq_reacIdx_dev,sizeof(int)*(nFromKeqStep*keqPadSize)),
     "hipMalloc(... fromKeq_reacIdx_dev ...)"
  );
  checkCudaError
  (
     hipMalloc((void**)&fromKeq_prodIdx_dev,sizeof(int)*(nFromKeqStep*keqPadSize)),
     "hipMalloc(... fromKeq_prodIdx_dev ...)"
  );
  checkCudaError
  (
     hipMalloc((void**)&fromKeq_stepIdx_dev,sizeof(int)*nFromKeqStep),
     "hipMalloc(... fromKeq_stepIdx_dev ...)"
  );
  checkCudaError
  (
     hipMalloc((void**)&fromKeq_nDelta_dev,sizeof(double)*nFromKeqStep),
     "hipMalloc(... fromKeq_nDelta_dev ...)"
  );
  checkCudaError
  (
     hipMalloc((void**)&fromKeq_stoich_fwd_dev,sizeof(double)*(nFromKeqStep*keqPadSize)),
     "hipMalloc(... fromKeq_stoich_fwd_dev ...)"
  );
  checkCudaError
  (
     hipMalloc((void**)&fromKeq_stoich_rev_dev,sizeof(double)*(nFromKeqStep*keqPadSize)),
     "hipMalloc(... fromKeq_stoich_rev_dev ...)"
  );
  hipMemcpy(fromKeq_reacIdx_dev,fromKeq_reacIdx_tmp,nFromKeqStep*keqPadSize*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(fromKeq_prodIdx_dev,fromKeq_prodIdx_tmp,nFromKeqStep*keqPadSize*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(fromKeq_stepIdx_dev,fromKeq_stepIdx_tmp,nFromKeqStep*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(fromKeq_nDelta_dev,fromKeq_nDelta_tmp,nFromKeqStep*sizeof(double),hipMemcpyHostToDevice);
  hipMemcpy(fromKeq_stoich_fwd_dev,fromKeq_stoich_fwd_tmp,nFromKeqStep*keqPadSize*sizeof(double),hipMemcpyHostToDevice);
  hipMemcpy(fromKeq_stoich_rev_dev,fromKeq_stoich_rev_tmp,nFromKeqStep*keqPadSize*sizeof(double),hipMemcpyHostToDevice);


  checkCudaError
  (
    hipMalloc((void**)&thirdBody_fwdStepIdx_dev,sizeof(int)*nThirdBodyRxn),
     "hipMalloc(... thirdBody_fwdStepIdx_dev ...)"
  );
  checkCudaError
  (
    hipMalloc((void**)&thirdBody_revStepIdx_dev,sizeof(int)*nThirdBodyRxn),
     "hipMalloc(... thirdBody_revStepIdx_dev ...)"
  );
  checkCudaError
  (
    hipMalloc((void**)&thirdBody_nEnhanced_dev,sizeof(int)*nThirdBodyRxn),
     "hipMalloc(... thirdBody_nEnhanced_dev ...)"
  );
  checkCudaError
  (
    hipMalloc((void**)&thirdBody_etbSpcIdx_dev,sizeof(int)*nThirdBodyRxn*maxThirdBodySpc),
     "hipMalloc(... thirdBody_etbSpcIdx_dev ...)"
  );
  checkCudaError
  (
    hipMalloc((void**)&thirdBody_etbSpcEff_dev,sizeof(double)*nThirdBodyRxn*maxThirdBodySpc),
     "hipMalloc(... thirdBody_etbSpcEff_dev ...)"
  );
  hipMemcpy(thirdBody_fwdStepIdx_dev,thirdBody_fwdStepIdx_tmp,nThirdBodyRxn*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(thirdBody_revStepIdx_dev,thirdBody_revStepIdx_tmp,nThirdBodyRxn*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(thirdBody_nEnhanced_dev,thirdBody_nEnhanced_tmp,nThirdBodyRxn*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(thirdBody_etbSpcIdx_dev,thirdBody_etbSpcIdx_tmp,nThirdBodyRxn*maxThirdBodySpc*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(thirdBody_etbSpcEff_dev,thirdBody_etbSpcEff_tmp,nThirdBodyRxn*maxThirdBodySpc*sizeof(double),hipMemcpyHostToDevice);


  checkCudaError
  (
    hipMalloc((void**)&falloff_falloffType_dev,sizeof(int)*nFalloffRxn),
     "hipMalloc(... falloff_falloffType_dev ...)"
  );
  checkCudaError
  (
    hipMalloc((void**)&falloff_fwdStepIdx_dev,sizeof(int)*nFalloffRxn),
     "hipMalloc(... falloff_fwdStepIdx_dev ...)"
  );
  checkCudaError
  (
    hipMalloc((void**)&falloff_revStepIdx_dev,sizeof(int)*nFalloffRxn),
     "hipMalloc(... falloff_revStepIdx_dev ...)"
  );
  checkCudaError
  (
    hipMalloc((void**)&falloff_nEnhanced_dev,sizeof(int)*nFalloffRxn),
     "hipMalloc(... falloff_nEnhanced_dev ...)"
  );
  checkCudaError
  (
    hipMalloc((void**)&falloff_falloffSpcIdx_dev,sizeof(int)*nFalloffRxn),
     "hipMalloc(... falloff_falloffSpcIdx_dev ...)"
  );
  checkCudaError
  (
    hipMalloc((void**)&falloff_etbSpcIdx_dev,sizeof(int)*nFalloffRxn*maxThirdBodySpc),
     "hipMalloc(... falloff_etbSpcIdx_dev ...)"
  );
  checkCudaError
  (
    hipMalloc((void**)&falloff_etbSpcEff_dev,sizeof(double)*nFalloffRxn*maxThirdBodySpc),
     "hipMalloc(... falloff_etbSpcEff_dev ...)"
  );
  checkCudaError
  (
    hipMalloc((void**)&falloff_param_dev,sizeof(double)*nFalloffRxn*maxFalloffParams),
     "hipMalloc(... falloff_param_dev ...)"
  );
  hipMemcpy(falloff_falloffType_dev,falloff_falloffType_tmp,nFalloffRxn*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(falloff_fwdStepIdx_dev,falloff_fwdStepIdx_tmp,nFalloffRxn*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(falloff_revStepIdx_dev,falloff_revStepIdx_tmp,nFalloffRxn*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(falloff_nEnhanced_dev,falloff_nEnhanced_tmp,nFalloffRxn*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(falloff_falloffSpcIdx_dev,falloff_falloffSpcIdx_tmp,nFalloffRxn*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(falloff_etbSpcIdx_dev,falloff_etbSpcIdx_tmp,nFalloffRxn*maxThirdBodySpc*sizeof(int),hipMemcpyHostToDevice);
  hipMemcpy(falloff_etbSpcEff_dev,falloff_etbSpcEff_tmp,nFalloffRxn*maxThirdBodySpc*sizeof(double),hipMemcpyHostToDevice);
  hipMemcpy(falloff_param_dev,falloff_param_tmp,nFalloffRxn*maxFalloffParams*sizeof(double),hipMemcpyHostToDevice);


  checkCudaError
  (
     hipMalloc((void**)&PLog_signAfact_dev,sizeof(double)*(nPLogInterpolationStep+nPLogInterpolationStepRev)*PLog_maxPressurePoints),
     "hipMalloc(... PLog_signAfact_dev ...)"
  );
  checkCudaError
  (
     hipMalloc((void**)&PLog_logAfact_dev,sizeof(double)*(nPLogInterpolationStep+nPLogInterpolationStepRev)*PLog_maxPressurePoints),
     "hipMalloc(... PLog_logAfact_dev ...)"
  );
  checkCudaError
  (
     hipMalloc((void**)&PLog_Tpow_dev,sizeof(double)*(nPLogInterpolationStep+nPLogInterpolationStepRev)*PLog_maxPressurePoints),
     "hipMalloc(... PLog_Tpow_dev ...)"
  );
  checkCudaError
  (
     hipMalloc((void**)&PLog_Tact_dev,sizeof(double)*(nPLogInterpolationStep+nPLogInterpolationStepRev)*PLog_maxPressurePoints),
     "hipMalloc(... PLog_Tact_dev ...)"
  );
  checkCudaError
  (
     hipMalloc((void**)&PLog_pressure_points_dev,sizeof(double)*(nPLogInterpolationStep+nPLogInterpolationStepRev)*PLog_maxPressurePoints),
     "hipMalloc(... PLog_pressure_points_dev ...)"
  );
  checkCudaError
  (
     hipMalloc((void**)&PLog_stepIdxs_dev,sizeof(int)*(nPLogInterpolationStep+nPLogInterpolationStepRev)),
     "hipMalloc(... PLog_stepIdxs_dev ...)"
  );
  hipMemcpy(PLog_signAfact_dev,&(plog_sign_Afact[0]),(nPLogInterpolationStep+nPLogInterpolationStepRev)*PLog_maxPressurePoints*sizeof(double),hipMemcpyHostToDevice);
  hipMemcpy(PLog_logAfact_dev,&(plog_logAfact[0]),(nPLogInterpolationStep+nPLogInterpolationStepRev)*PLog_maxPressurePoints*sizeof(double),hipMemcpyHostToDevice);
  hipMemcpy(PLog_Tpow_dev,&(plog_Tpow[0]),(nPLogInterpolationStep+nPLogInterpolationStepRev)*PLog_maxPressurePoints*sizeof(double),hipMemcpyHostToDevice);
  hipMemcpy(PLog_Tact_dev,&(plog_Tact[0]),(nPLogInterpolationStep+nPLogInterpolationStepRev)*PLog_maxPressurePoints*sizeof(double),hipMemcpyHostToDevice);
  hipMemcpy(PLog_pressure_points_dev,&(plog_pressure_points[0]),(nPLogInterpolationStep+nPLogInterpolationStepRev)*PLog_maxPressurePoints*sizeof(double),hipMemcpyHostToDevice);
  hipMemcpy(PLog_stepIdxs_dev,&(plog_stepIdxs[0]),(nPLogInterpolationStep+nPLogInterpolationStepRev)*sizeof(int),hipMemcpyHostToDevice);


  hipStreamCreate(&arrhStream);
  hipStreamCreate(&kEqStream);
  hipStreamCreate(&thirdBodyStream);
  hipStreamCreate(&falloffStream);
  hipStreamCreate(&PLogStream);

#ifdef ZERORK_FULL_DEBUG
  hipDeviceSynchronize();
  checkCudaError(hipGetLastError(),"rate_const_cuda::setGpuParams()");
#endif
}


void rate_const_cuda::updateK_CUDA_mr
(
    const int nReactors, const double * T_dev, const double * C_dev, double * K_dev
)
{
  rate_const_updateArrheniusStep_CUDA_mr(nReactors,nStep,K_dev,logAfact_dev,Tpow_dev,Tact_dev,T_dev,arrhStream);
  rate_const_concentrationSum(nReactors, nSpc, C_dev, Csum_dev, PLogStream);

  if(nFromKeqStep > 0) {
    hipMemsetAsync(&Gibbs_RT_dev[nSpc*nReactors],0,sizeof(double)*nReactors,kEqStream);//N.B. this is zero-ing the padded space
    static_cast<nasa_poly_group_cuda*>(thermoPtr)->getG_RT_CUDA_mr(nReactors,T_dev,Gibbs_RT_dev,kEqStream);
  }

  rate_const_updatePLogInterpolationStep_CUDA_mr(nReactors,
          nPLogInterpolationStep+nPLogInterpolationStepRev, PLog_maxPressurePoints,
          PLog_stepIdxs_dev, PLog_pressure_points_dev, PLog_signAfact_dev,
          PLog_logAfact_dev, PLog_Tpow_dev, PLog_Tact_dev, Csum_dev, T_dev,
          K_dev, PLogStream);

  //N.B. Keq is done in default stream which is synchronous across all streams
  if(!use_non_integer_network_) {
    rate_const_updateFromKeqStep_CUDA_mr(nReactors, nSpc, nStep, nFromKeqStep, keqPadSize,
            fromKeq_reacIdx_dev, fromKeq_prodIdx_dev, fromKeq_stepIdx_dev,
            fromKeq_nDelta_dev, Gibbs_RT_dev, K_dev, T_dev, (hipStream_t) 0);
  } else {
    rate_const_updateFromKeqStepNI_CUDA_mr(nReactors, nSpc, nStep, nFromKeqStep, keqPadSize,
            fromKeq_reacIdx_dev, fromKeq_prodIdx_dev, fromKeq_stepIdx_dev,
            fromKeq_nDelta_dev, fromKeq_stoich_fwd_dev, fromKeq_stoich_rev_dev,
            Gibbs_RT_dev, K_dev, T_dev, (hipStream_t) 0);
  }
  if(nFromKeqStep == 0) {
    hipStreamSynchronize(arrhStream);
    hipStreamSynchronize(PLogStream);
  }

  rate_const_updateFalloff_CUDA_mr(nReactors, nSpc, nStep, nFalloffRxn, maxThirdBodySpc,
         maxFalloffParams, falloff_falloffType_dev, falloff_nEnhanced_dev,
         falloff_falloffSpcIdx_dev,
         falloff_etbSpcIdx_dev, falloff_fwdStepIdx_dev, falloff_revStepIdx_dev,
         falloff_etbSpcEff_dev, falloff_param_dev, logAfact_dev, Tpow_dev, Tact_dev,
         C_dev, Csum_dev, T_dev, K_dev, falloffStream);

  rate_const_updateThirdBody_CUDA_mr(nReactors, nSpc, nStep, nThirdBodyRxn, maxThirdBodySpc,
         thirdBody_nEnhanced_dev, thirdBody_etbSpcIdx_dev, thirdBody_fwdStepIdx_dev,
         thirdBody_revStepIdx_dev, thirdBody_etbSpcEff_dev, C_dev, Csum_dev, K_dev,
         thirdBodyStream);
}

} // namespace zerork

