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

  void calcRatesFromTC_CUDA(const double T, const double C[],
           double netOut[], double createOut[], double destroyOut[],
     double stepOut[]);
  void calcRatesFromTC_CUDA_mr(const int nReactors, const double T[], const double C[],
           double netOut[], double createOut[], double destroyOut[],
     double stepOut[], bool transposeInOut = false, bool outputAll = false);
  void calcRatesFromTC_CUDA_mr_dev(const int nReactors, const double T_dev[], const double C_dev[],
           double *netOut_dev, double *createOut_dev, double *destroyOut_dev,
     double *stepOut_dev);
  int reorderScatterAdd_by_ntuple(const int ntuple, const int nOps, const int srcSize,
                     const int destSize, int srcId[], int destId[]);

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
  double *stepOut_dev;
  double *createOut_dev;
  double *destroyOut_dev;
  double *netOut_dev;
  double *Tmulti_dev;

  int scatterAddCount;
  cudaStream_t * scatterAddStreams;
  //scatter add starting indicies
  int sa64c,sa32c,sa16c,sa8c,sa4c,sa2c,sa1c;
  int sa64d,sa32d,sa16d,sa8d,sa4d,sa2d,sa1d;
  int maxOpsDestroy;
  int maxOpsCreate;

  double gpuKTime,gpuStepTime,gpuProdTime,gpuNetTime,gpuTxTime;
};

} // namespace zerork

#endif
