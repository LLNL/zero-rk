#ifndef ZERORK_PERF_NET_CUDA_H
#define ZERORK_PERF_NET_CUDA_H

#include "info_net.h"
#include "rate_const_cuda.h"
#include "perf_net.h"
#include "external_funcs.h"
#include <cuda_runtime.h>

namespace zerork {


class perf_net_cuda : public perf_net
{
 public: 
  perf_net_cuda(info_net &netobj, rate_const &Kobj);
  virtual ~perf_net_cuda();

  void calcRatesFromTC_CUDA_mr_dev(const int nReactors, const double T_dev[],
		                   const double C_dev[], const double stepLimiter_dev[],
				   double *netOut_dev, double *createOut_dev, double *destroyOut_dev,
                                   double *stepOut_dev);
  int reorderScatterAdd_by_ntuple(const int ntuple, const int nOps, const int srcSize,
                     const int destSize, int srcId[], int destId[], double* srcMult=nullptr);

 private:
  int *reactantSpcIdxListUnwrapped;
  int *reactantSpcIdxListUnwrapped_dev;
  int *reactantSpcIdxList_dev;
  int *reactantStepIdxList_dev;
  int *productSpcIdxList_dev;
  int *productStepIdxList_dev;
  int *nOpsDestroy_dev;
  int *nOpsCreate_dev;
  double *C_dev;

  double* rop_concentration_powers_dev;
  int *niReactantSpcIdxList_dev;
  int *niReactantStepIdxList_dev;
  double *niReactantStoichNumList_dev;
  int *niProductSpcIdxList_dev;
  int* niProductStepIdxList_dev;
  double *niProductStoichNumList_dev;
  int* nOpsDestroy_ni_dev;
  int* nOpsCreate_ni_dev;


  int scatterAddCount;
  cudaStream_t * scatterAddStreams;
  //scatter add starting indicies
  int sa64c,sa32c,sa16c,sa8c,sa4c,sa2c,sa1c;
  int sa64d,sa32d,sa16d,sa8d,sa4d,sa2d,sa1d;
  int maxOpsDestroy;
  int maxOpsCreate;

  int sa64c_ni,sa32c_ni,sa16c_ni,sa8c_ni,sa4c_ni,sa2c_ni,sa1c_ni;
  int sa64d_ni,sa32d_ni,sa16d_ni,sa8d_ni,sa4d_ni,sa2d_ni,sa1d_ni;
  int maxOpsDestroy_ni;
  int maxOpsCreate_ni;

  double gpuKTime,gpuStepTime,gpuProdTime,gpuNetTime,gpuTxTime;
};

} // namespace zerork

#endif
