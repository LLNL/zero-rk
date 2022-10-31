#include <cassert>
#include "perf_net_cuda.h"
#include <algorithm>

#include "zerork_cuda_defs.h"
#include "perf_net_kernels.h"
#include "misc_kernels.h"

#include <sys/time.h>

//#define ZERORK_CUDA_EVENTS
//#define ZERORK_SCATTER_ADD_FUSED

namespace zerork {

static double getHighResolutionTime(void)
{
    struct timeval tod; 

    gettimeofday(&tod, NULL);
    double time_seconds = (double) tod.tv_sec + ((double) tod.tv_usec / 1000000.0);
    return time_seconds;
}


perf_net_cuda::perf_net_cuda(info_net &netobj, rate_const &Kobj)
  : perf_net(netobj, Kobj)
{
  int j,k;
  int rCtrU;
  int nReactorsMax = static_cast<rate_const_cuda*>(rateConstPtr)->nReactorsMax();

  reactantSpcIdxListUnwrapped = new int[nStep*maxReactants];

  rCtrU = 0; // reactant species counter (for unwrapped list)
  for(j=0; j<nStep; j++)
  {
    for(k=0; k<netobj.getOrderOfStep(j); k++) 
    {    
      reactantSpcIdxListUnwrapped[rCtrU]  = netobj.getSpecIdxOfStepReactant(j,k);
      ++rCtrU;
    }    
    for(k=netobj.getOrderOfStep(j); k < maxReactants; k++) 
    {    
      reactantSpcIdxListUnwrapped[rCtrU]  = nSpc;
      ++rCtrU;
    }
  }
  assert(rCtrU == nStep*maxReactants);

  //Reorder scatter add operations to reduce atomic races in production rate calcs.
  sa64c=reorderScatterAdd_by_ntuple(128, totProd, nStep, nSpc,
                                       productStepIdxList,
                                       productSpcIdxList);
  maxOpsCreate = std::max(sa64c/128,0);

  sa64d=reorderScatterAdd_by_ntuple(128, totReac, nStep, nSpc,
                                       reactantStepIdxList,
                                       reactantSpcIdxList);
  maxOpsDestroy = std::max(sa64d/128,0);

  sa32c=reorderScatterAdd_by_ntuple(64, totProd-sa64c, nStep, nSpc,
                                       &productStepIdxList[sa64c],
                                       &productSpcIdxList[sa64c]);
  maxOpsCreate = std::max(sa32c/64,maxOpsCreate);
  sa32c+=sa64c;

  sa32d=reorderScatterAdd_by_ntuple(64, totReac-sa64d, nStep, nSpc,
                                       &reactantStepIdxList[sa64d],
                                       &reactantSpcIdxList[sa64d]);
  maxOpsDestroy = std::max(sa32d/64,maxOpsDestroy);
  sa32d+=sa64d;

  sa16c=reorderScatterAdd_by_ntuple(32, totProd-sa32c, nStep, nSpc,
                                       &productStepIdxList[sa32c],
                                       &productSpcIdxList[sa32c]);
  maxOpsCreate = std::max(sa16c/32,maxOpsCreate);
  sa16c+=sa32c;

  sa16d=reorderScatterAdd_by_ntuple(32, totReac-sa32d, nStep, nSpc,
                                       &reactantStepIdxList[sa32d],
                                       &reactantSpcIdxList[sa32d]);
  maxOpsDestroy = std::max(sa16d/32,maxOpsDestroy);
  sa16d+=sa32d;

  sa8c=reorderScatterAdd_by_ntuple(16, totProd-sa16c, nStep, nSpc,
                                       &productStepIdxList[sa16c],
                                       &productSpcIdxList[sa16c]);
  maxOpsCreate = std::max(sa8c/16,maxOpsCreate);
  sa8c+=sa16c;

  sa8d=reorderScatterAdd_by_ntuple(16, totReac-sa16d, nStep, nSpc,
                                       &reactantStepIdxList[sa16d],
                                       &reactantSpcIdxList[sa16d]);
  maxOpsDestroy = std::max(sa8d/16,maxOpsDestroy);
  sa8d+=sa16d;

  sa4c=reorderScatterAdd_by_ntuple(8, totProd-sa8c, nStep, nSpc,
                                       &productStepIdxList[sa8c],
                                       &productSpcIdxList[sa8c]);
  maxOpsCreate = std::max(sa4c/8,maxOpsCreate);
  sa4c+=sa8c;

  sa4d=reorderScatterAdd_by_ntuple(8, totReac-sa8d, nStep, nSpc,
                                       &reactantStepIdxList[sa8d],
                                       &reactantSpcIdxList[sa8d]);
  maxOpsDestroy = std::max(sa4d/8,maxOpsDestroy);
  sa4d+=sa8d;

  sa2c=reorderScatterAdd_by_ntuple(4, totProd-sa4c, nStep, nSpc,
                                       &productStepIdxList[sa4c],
                                       &productSpcIdxList[sa4c]);
  maxOpsCreate = std::max(sa2c/4,maxOpsCreate);
  sa2c+=sa4c;

  sa2d=reorderScatterAdd_by_ntuple(4, totReac-sa4d, nStep, nSpc,
                                       &reactantStepIdxList[sa4d],
                                       &reactantSpcIdxList[sa4d]);
  maxOpsDestroy = std::max(sa2d/4,maxOpsDestroy);
  sa2d+=sa4d;

  sa1c=reorderScatterAdd_by_ntuple(2, totProd-sa2c, nStep, nSpc,
                                       &productStepIdxList[sa2c],
                                       &productSpcIdxList[sa2c]);
  maxOpsCreate = std::max(sa1c/2,maxOpsCreate);
  sa1c+=sa2c;

  sa1d=reorderScatterAdd_by_ntuple(2, totReac-sa2d, nStep, nSpc,
                                       &reactantStepIdxList[sa2d],
                                       &reactantSpcIdxList[sa2d]);
  maxOpsDestroy = std::max(sa1d/2,maxOpsDestroy);
  sa1d+=sa2d;
  maxOpsCreate  = std::max(totProd-sa1c,maxOpsCreate);
  maxOpsDestroy = std::max(totReac-sa1d,maxOpsDestroy);

  int* nOpsDestroy = new int[9];
  nOpsDestroy[8] = 0;
  nOpsDestroy[7] = sa64d;
  nOpsDestroy[6] = sa32d;
  nOpsDestroy[5] = sa16d;
  nOpsDestroy[4] = sa8d;
  nOpsDestroy[3] = sa4d;
  nOpsDestroy[2] = sa2d;
  nOpsDestroy[1] = sa1d;
  nOpsDestroy[0] = totReac;

  int* nOpsCreate = new int[9];
  nOpsCreate[8] = 0;
  nOpsCreate[7] = sa64c;
  nOpsCreate[6] = sa32c;
  nOpsCreate[5] = sa16c;
  nOpsCreate[4] = sa8c;
  nOpsCreate[3] = sa4c;
  nOpsCreate[2] = sa2c;
  nOpsCreate[1] = sa1c;
  nOpsCreate[0] = totProd;

  //8 create and 8 destroy kernels  2^[0-7]
  scatterAddCount = 16;
  scatterAddStreams = new cudaStream_t[scatterAddCount];
  for(j = 0; j < scatterAddCount; ++j)
  {
      cudaStreamCreate( &scatterAddStreams[j] );
  }

  checkCudaError
  (
     cudaMalloc((void**)&nOpsDestroy_dev,sizeof(int)*9),
     "cudaMalloc(... nOpsDestroy_dev ...)"
  );
  cudaMemcpy(nOpsDestroy_dev, nOpsDestroy,sizeof(int)*9,cudaMemcpyHostToDevice);
  delete[] nOpsDestroy;

  checkCudaError
  (
     cudaMalloc((void**)&nOpsCreate_dev,sizeof(int)*9),
     "cudaMalloc(... nOpsCreate_dev ...)"
  );
  cudaMemcpy(nOpsCreate_dev, nOpsCreate,sizeof(int)*9,cudaMemcpyHostToDevice);
  delete[] nOpsCreate;

  checkCudaError
  (
     cudaMalloc((void**)&reactantSpcIdxListUnwrapped_dev,sizeof(int)*rCtrU),
     "cudaMalloc(... reactantSpcIdxListUnwrapped_dev ...)"
  );
  cudaMemcpy(reactantSpcIdxListUnwrapped_dev,reactantSpcIdxListUnwrapped,sizeof(int)*rCtrU,cudaMemcpyHostToDevice);

  checkCudaError
  (
     cudaMalloc((void**)&reactantSpcIdxList_dev,sizeof(int)*totReac),
     "cudaMalloc(... reactantSpcIdxList_dev ...)"
  );
  cudaMemcpy(reactantSpcIdxList_dev,reactantSpcIdxList,sizeof(int)*totReac,cudaMemcpyHostToDevice);
  checkCudaError
  (
     cudaMalloc((void**)&reactantStepIdxList_dev,sizeof(int)*totReac),
     "cudaMalloc(... reactantStepIdxList_dev ...)"
  );
  cudaMemcpy(reactantStepIdxList_dev,reactantStepIdxList,sizeof(int)*totReac,cudaMemcpyHostToDevice);


  checkCudaError
  (
     cudaMalloc((void**)&productSpcIdxList_dev,sizeof(int)*totProd),
     "cudaMalloc(... productSpcIdxList_dev ...)"
  );
  cudaMemcpy(productSpcIdxList_dev,productSpcIdxList,sizeof(int)*totProd,cudaMemcpyHostToDevice);
  checkCudaError
  (
     cudaMalloc((void**)&productStepIdxList_dev,sizeof(int)*totProd),
     "cudaMalloc(... productStepIdxList_dev ...)"
  );
  cudaMemcpy(productStepIdxList_dev,productStepIdxList,sizeof(int)*totProd,cudaMemcpyHostToDevice);

  checkCudaError
  (
     cudaMalloc((void**)&stepOut_dev,sizeof(double)*nStep*nReactorsMax),
     "cudaMalloc(... stepOut_dev ...)"
  );

  checkCudaError
  (
     cudaMalloc((void**)&createOut_dev,sizeof(double)*nSpc*nReactorsMax),
     "cudaMalloc(... createOut_dev ...)"
  );
  checkCudaError
  (
     cudaMalloc((void**)&destroyOut_dev,sizeof(double)*nSpc*nReactorsMax),
     "cudaMalloc(... destroyOut_dev ...)"
  );
  checkCudaError
  (
     cudaMalloc((void**)&netOut_dev,sizeof(double)*nSpc*nReactorsMax),
     "cudaMalloc(... netOut_dev ...)"
  );

  checkCudaError
  (
     cudaMalloc((void**)&C_dev,sizeof(double)*(nSpc+1)*nReactorsMax),
     "cudaMalloc(... C_dev ...)"
  );
  //Set extra concentration to 1.0 to keep scatter mult kernel in sync
  setDoubleArrayVal(C_dev,1.0,nSpc*nReactorsMax,(nSpc+1)*nReactorsMax);

  checkCudaError
  (
     cudaMalloc((void**)&Tmulti_dev,sizeof(double)*nReactorsMax),
     "cudaMalloc(... Tmulti_dev ...)"
  );

  //Initialize gpu timing data
  gpuKTime = gpuStepTime = gpuProdTime = gpuNetTime = gpuTxTime = 0.0;
}

perf_net_cuda::~perf_net_cuda()
{
#ifdef ZERORK_CUDA_EVENTS
  double cpuTotTime = cpuKTime + cpuStepTime + cpuProdTime + cpuNetTime;
  double gpuTotTime = gpuKTime + gpuStepTime + gpuProdTime + gpuNetTime + gpuTxTime;
  printf("#Timing data: KTime         StepTime      ProdTime      NetTime       TransferTime  TotalTime    \n");
  printf("#  CPU:       %11.7e %11.7e %11.7e %11.7e %11.7e %11.7e\n",
                   cpuKTime,cpuStepTime,cpuProdTime,cpuNetTime,0.0,cpuTotTime);
  printf("#  GPU:       %11.7e %11.7e %11.7e %11.7e %11.7e %11.7e\n",
                   gpuKTime,gpuStepTime,gpuProdTime,gpuNetTime,gpuTxTime,gpuTotTime);
#endif

  delete [] reactantSpcIdxListUnwrapped;
  cudaFree(reactantSpcIdxListUnwrapped_dev);
  cudaFree(nOpsDestroy_dev);
  cudaFree(nOpsCreate_dev);
  cudaFree(stepOut_dev);
  cudaFree(createOut_dev);
  cudaFree(destroyOut_dev);
  cudaFree(netOut_dev);
  cudaFree(C_dev);
  cudaFree(productSpcIdxList_dev);
  cudaFree(productStepIdxList_dev);
  cudaFree(reactantSpcIdxList_dev);
  cudaFree(reactantStepIdxList_dev);
  cudaFree(Tmulti_dev);
  for(int j = 0; j < scatterAddCount; ++j)
  {
      cudaStreamDestroy( scatterAddStreams[j] );
  }
  delete [] scatterAddStreams;
}

void perf_net_cuda::calcRatesFromTC_CUDA(const double T, const double C[],
				   double netOut[], double createOut[],
				   double destroyOut[], double stepOut[])
{
    int cpySpcSize = nSpc*sizeof(double);
    int cpyStepSize = nStep*sizeof(double);
    float cudaTime_ms;
#ifdef ZERORK_CUDA_EVENTS
    cudaEvent_t startTxC,endTxC,startTxStep,endTxStep,startK,endK,startStep,endStep,startProd,endProd,startNet,endNet;
    cudaEvent_t startTxProd,endTxProd,startTxNet,endTxNet;
    cudaEventCreate(&startTxC);
    cudaEventCreate(&endTxC);
    cudaEventCreate(&startTxStep);
    cudaEventCreate(&endTxStep);
    cudaEventCreate(&startK);
    cudaEventCreate(&endK);
    cudaEventCreate(&startStep);
    cudaEventCreate(&endStep);
    cudaEventCreate(&startProd);
    cudaEventCreate(&endProd);
    cudaEventCreate(&startNet);
    cudaEventCreate(&endNet);
    cudaEventCreate(&startTxProd);
    cudaEventCreate(&endTxProd);
    cudaEventCreate(&startTxNet);
    cudaEventCreate(&endTxNet);
#endif

    // reaction-species link lists are kept on device with work arrays for 
    // rxn, creation, destruction, and net species rates
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(startTxC);
#endif
    cudaMemcpy(C_dev,C,cpySpcSize,cudaMemcpyHostToDevice);
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endTxC);
#endif

    setDoubleArrayVal(C_dev,1.0,nSpc,nSpc+1);

#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(startK);
#endif
    static_cast<rate_const_cuda*>(rateConstPtr)->updateK_CUDA(T,C,C_dev,stepOut_dev);
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endK);
#endif

#define ZERORK_GPU_SCATTER_MULT
#ifdef ZERORK_GPU_SCATTER_MULT
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(startStep);
#endif
    //Scatter multiplication kernel
    perf_net_cuda_rxn_conc_mult(nStep, maxReactants, reactantSpcIdxListUnwrapped_dev, C_dev, stepOut_dev);
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endStep);
#endif

#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(startTxStep);
#endif
    cudaMemcpyAsync(stepOut,stepOut_dev,cpyStepSize,cudaMemcpyDeviceToHost);
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endTxStep);
#endif
#else
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(startTxStep);
#endif
    cudaMemcpy(stepOut,stepOut_dev,cpyStepSize,cudaMemcpyDeviceToHost);
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endTxStep);
#endif

#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(startStep);
#endif
    for(int j=0; j<totReac; ++j)
      {stepOut[reactantStepIdxList[j]]*=C[reactantSpcIdxList[j]];}
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endStep);
#endif
#endif
#ifdef ZERORK_CUDA_EVENTS
    cudaEventSynchronize(endTxC);
    cudaEventSynchronize(endK);
    cudaEventSynchronize(endStep);

    cudaEventElapsedTime(&cudaTime_ms, startTxC, endTxC);
    gpuTxTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startK, endK);
    gpuKTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startStep, endStep);
    gpuStepTime += cudaTime_ms/1000.0;
#endif


#define GPU_DO_ALL
#ifdef GPU_DO_ALL
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(startProd);
#endif
    cudaMemset(createOut_dev,0,nSpc*sizeof(double));
    cudaMemset(destroyOut_dev,0,nSpc*sizeof(double));

    perf_net_scatterAdd_gpu_atomic_global_128op
        (sa64d,
         &reactantStepIdxList_dev[0],
         &reactantSpcIdxList_dev[0],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         destroyOut_dev, scatterAddStreams[0]);
      perf_net_scatterAdd_gpu_atomic_global_64op
        (sa32d-sa64d,
         &reactantStepIdxList_dev[sa64d],
         &reactantSpcIdxList_dev[sa64d],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         destroyOut_dev, scatterAddStreams[1]);
      perf_net_scatterAdd_gpu_atomic_global_32op
        (sa16d-sa32d,
         &reactantStepIdxList_dev[sa32d],
         &reactantSpcIdxList_dev[sa32d],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         destroyOut_dev, scatterAddStreams[2]);
      perf_net_scatterAdd_gpu_atomic_global_16op
        (sa8d-sa16d,
         &reactantStepIdxList_dev[sa16d],
         &reactantSpcIdxList_dev[sa16d],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         destroyOut_dev, scatterAddStreams[3]);
      perf_net_scatterAdd_gpu_atomic_global_8op
        (sa4d-sa8d,
         &reactantStepIdxList_dev[sa8d],
         &reactantSpcIdxList_dev[sa8d],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         destroyOut_dev, scatterAddStreams[4]);
       perf_net_scatterAdd_gpu_atomic_global_4op
        (sa2d-sa4d,
         &reactantStepIdxList_dev[sa4d],
         &reactantSpcIdxList_dev[sa4d],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         destroyOut_dev, scatterAddStreams[5]);
      perf_net_scatterAdd_gpu_atomic_global_2op
        (sa1d-sa2d,
         &reactantStepIdxList_dev[sa2d],
         &reactantSpcIdxList_dev[sa2d],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         destroyOut_dev, scatterAddStreams[6]);
      perf_net_scatterAdd_gpu_atomic_global
        (totReac-sa1d,
         &reactantStepIdxList_dev[sa1d],
         &reactantSpcIdxList_dev[sa1d],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         destroyOut_dev, scatterAddStreams[7]);

      perf_net_scatterAdd_gpu_atomic_global_128op
        (sa64c,
         &productStepIdxList_dev[0],
         &productSpcIdxList_dev[0],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         createOut_dev, scatterAddStreams[8]);
      perf_net_scatterAdd_gpu_atomic_global_64op
        (sa32c-sa64c,
         &productStepIdxList_dev[sa64c],
         &productSpcIdxList_dev[sa64c],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         createOut_dev, scatterAddStreams[9]);
      perf_net_scatterAdd_gpu_atomic_global_32op
        (sa16c-sa32c,
         &productStepIdxList_dev[sa32c],
         &productSpcIdxList_dev[sa32c],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         createOut_dev, scatterAddStreams[10]);
      perf_net_scatterAdd_gpu_atomic_global_16op
        (sa8c-sa16c,
         &productStepIdxList_dev[sa16c],
         &productSpcIdxList_dev[sa16c],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         createOut_dev, scatterAddStreams[11]);
      perf_net_scatterAdd_gpu_atomic_global_8op
        (sa4c-sa8c,
         &productStepIdxList_dev[sa8c],
         &productSpcIdxList_dev[sa8c],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         createOut_dev, scatterAddStreams[12]);
       perf_net_scatterAdd_gpu_atomic_global_4op
        (sa2c-sa4c,
         &productStepIdxList_dev[sa4c],
         &productSpcIdxList_dev[sa4c],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         createOut_dev, scatterAddStreams[13]);
      perf_net_scatterAdd_gpu_atomic_global_2op
        (sa1c-sa2c,
         &productStepIdxList_dev[sa2c],
         &productSpcIdxList_dev[sa2c],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         createOut_dev, scatterAddStreams[14]);
      perf_net_scatterAdd_gpu_atomic_global
        (totProd-sa1c,
         &productStepIdxList_dev[sa1c],
         &productSpcIdxList_dev[sa1c],
         1,
         nStep,
         stepOut_dev,
         nSpc,
         createOut_dev, scatterAddStreams[15]);

#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endProd);
    cudaEventSynchronize(endProd);
#else
    for(int j = 0; j <16; ++j) {
      cudaStreamSynchronize(scatterAddStreams[j]);
    }
#endif

#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(startTxProd);
#endif
    cudaMemcpyAsync(createOut,createOut_dev,cpySpcSize,cudaMemcpyDeviceToHost);
    cudaMemcpyAsync(destroyOut,destroyOut_dev,cpySpcSize,cudaMemcpyDeviceToHost);
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endTxProd);
#endif

#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(startNet);
#endif
    perf_net_cuda_net_rates(nSpc, createOut_dev, destroyOut_dev, netOut_dev);
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endNet);
#endif

#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(startTxNet);
#endif
    cudaMemcpy(netOut,netOut_dev,cpySpcSize,cudaMemcpyDeviceToHost);
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endTxNet);
#endif

#ifdef ZERORK_CUDA_EVENTS
    cudaEventSynchronize(endTxProd);
    cudaEventSynchronize(endTxNet);


    cudaEventElapsedTime(&cudaTime_ms, startProd, endProd);
    gpuProdTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startTxProd, endTxProd);
    gpuTxTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startNet, endNet);
    gpuNetTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startTxNet, endTxNet);
    gpuTxTime += cudaTime_ms/1000.0;
#endif
#else //GPU_DO_ALL

    double startTime = getHighResolutionTime();
    memset(createOut,0,nSpc*sizeof(double));
    memset(destroyOut,0,nSpc*sizeof(double));

#ifdef ZERORK_CUDA_EVENTS
    cudaEventSynchronize(endTxStep);
    cudaEventElapsedTime(&cudaTime_ms, startTxStep, endTxStep);
    gpuTxTime += cudaTime_ms/1000.0;
#endif
    // compute the species destruction rate by adding each steps rate of progress
    // to the sum for each reactant species found
    for(int j=0; j<totReac; ++j)
      {destroyOut[reactantSpcIdxList[j]]+=stepOut[reactantStepIdxList[j]];}

    // compute the species creation rate by adding each steps rate of progress
    // to the sum for each product species found 
    for(int j=0; j<totProd; ++j)
      {createOut[productSpcIdxList[j]]+=stepOut[productStepIdxList[j]];}
    gpuProdTime += getHighResolutionTime() - startTime;

    startTime = getHighResolutionTime();
    // compute the net species production rate = create - destroy
    int loopLim = nSpc;
    for(int j=0; j<loopLim; ++j )
      {netOut[j]=createOut[j]-destroyOut[j];}
    gpuNetTime += getHighResolutionTime() - startTime;
#endif
#ifdef ZERORK_CUDA_EVENTS
    cudaEventDestroy(startTxC);
    cudaEventDestroy(endTxC);
    cudaEventDestroy(startTxStep);
    cudaEventDestroy(endTxStep);
    cudaEventDestroy(startK);
    cudaEventDestroy(endK);
    cudaEventDestroy(startStep);
    cudaEventDestroy(endStep);
    cudaEventDestroy(startProd);
    cudaEventDestroy(endProd);
    cudaEventDestroy(startNet);
    cudaEventDestroy(endNet);
    cudaEventDestroy(startTxProd);
    cudaEventDestroy(endTxProd);
    cudaEventDestroy(startTxNet);
    cudaEventDestroy(endTxNet);
#endif
}


void perf_net_cuda::calcRatesFromTC_CUDA_mr(const int nReactors, const double T[], const double C[],
				   double netOut[], double createOut[],
				   double destroyOut[], double stepOut[],bool transposeInOut, bool outputAll)
{
    int j,k;//,spcOffset,stepOffset;
    int cpySpcSize = nReactors*nSpc*sizeof(double);
    int cpyStepSize = nReactors*nStep*sizeof(double);
    float cudaTime_ms;
    double *netOut_tp,*stepOut_tp,*createOut_tp,*destroyOut_tp;

    //Reset dummy values in case nReactors is different from last call.
    setDoubleArrayVal(C_dev,1.0,nSpc*nReactors,(nSpc+1)*nReactors);

    if(transposeInOut)
    {
        netOut_tp = (double*)malloc(sizeof(double)*nSpc*nReactors);
        if(outputAll)
        {
            createOut_tp = (double*)malloc(sizeof(double)*nSpc*nReactors);
            destroyOut_tp = (double*)malloc(sizeof(double)*nSpc*nReactors);
            stepOut_tp = (double*)malloc(sizeof(double)*nStep*nReactors);
        }
    }
    else
    {
        netOut_tp = netOut;
        createOut_tp = createOut;
        destroyOut_tp = destroyOut;
        stepOut_tp = stepOut;
    }

    cudaEvent_t startTxC,endTxC,startTxStep,endTxStep,startK,endK,startStep,endStep,
                startTxProd,endTxProd,startProd,endProd,startTxNet,endTxNet,startNet,endNet;
    cudaEventCreate(&startTxC);
    cudaEventCreate(&endTxC);
    cudaEventCreate(&startTxStep);
    cudaEventCreate(&endTxStep);
    cudaEventCreate(&startK);
    cudaEventCreate(&endK);
    cudaEventCreate(&startStep);
    cudaEventCreate(&endStep);
    cudaEventCreate(&startTxProd);
    cudaEventCreate(&endTxProd);
    cudaEventCreate(&startProd);
    cudaEventCreate(&endProd);
    cudaEventCreate(&startTxNet);
    cudaEventCreate(&endTxNet);
    cudaEventCreate(&startNet);
    cudaEventCreate(&endNet);

    // transpose C
    if(transposeInOut)
    {
        double *C_tp = new double[nSpc*nReactors];
        for(j=0;j < nReactors; ++j)
        {
            for(k=0;k<nSpc;++k)
            {
                int origIdx = (j*nSpc+k);
                int transIdx = (k*nReactors+j);
                C_tp[transIdx] = C[origIdx];
            }
        }
        // reaction-species link lists are kept on device with work arrays for 
        // rxn, creation, destruction, and net species rates
        cudaEventRecord(startTxC);
        cudaMemcpy(C_dev,C_tp,cpySpcSize,cudaMemcpyHostToDevice);
        cudaEventRecord(endTxC);
        cudaDeviceSynchronize();
        delete [] C_tp;
    }
    else
    {
        cudaEventRecord(startTxC);
        cudaMemcpy(C_dev,C,cpySpcSize,cudaMemcpyHostToDevice);
        cudaEventRecord(endTxC);
    }

    //Do here instead of in updateK to keep from having separate updateK for tp_dev
    cudaMemcpy(Tmulti_dev,T,sizeof(double)*nReactors,cudaMemcpyHostToDevice);

    cudaEventRecord(startK);
    static_cast<rate_const_cuda*>(rateConstPtr)->updateK_CUDA_mr(nReactors,Tmulti_dev,C_dev,stepOut_dev);
    cudaEventRecord(endK);

    cudaEventRecord(startStep);
    perf_net_cuda_rxn_conc_mult_mr(nReactors, nSpc, nStep,
        maxReactants, reactantSpcIdxListUnwrapped_dev, C_dev, stepOut_dev);
    cudaEventRecord(endStep);

    cudaEventRecord(startTxStep);
    if(outputAll)
    {
        cudaMemcpyAsync(stepOut_tp,stepOut_dev,cpyStepSize,cudaMemcpyDeviceToHost);
    }
    cudaEventRecord(endTxStep);

    cudaEventSynchronize(endTxC);
    cudaEventSynchronize(endK);
    cudaEventSynchronize(endStep);

    cudaEventElapsedTime(&cudaTime_ms, startTxC, endTxC);
    gpuTxTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startK, endK);
    gpuKTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startStep, endStep);
    gpuStepTime += cudaTime_ms/1000.0;

    cudaEventRecord(startProd);
    cudaMemset(createOut_dev,0,cpySpcSize);
    cudaMemset(destroyOut_dev,0,cpySpcSize);


    perf_net_scatterAdd_gpu_atomic_global_128op
            (sa64d,
             &reactantStepIdxList_dev[0],
             &reactantSpcIdxList_dev[0],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             destroyOut_dev, scatterAddStreams[0]);
          perf_net_scatterAdd_gpu_atomic_global_64op
            (sa32d-sa64d,
             &reactantStepIdxList_dev[sa64d],
             &reactantSpcIdxList_dev[sa64d],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             destroyOut_dev, scatterAddStreams[1]);
          perf_net_scatterAdd_gpu_atomic_global_32op
            (sa16d-sa32d,
             &reactantStepIdxList_dev[sa32d],
             &reactantSpcIdxList_dev[sa32d],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             destroyOut_dev, scatterAddStreams[2]);
          perf_net_scatterAdd_gpu_atomic_global_16op
            (sa8d-sa16d,
             &reactantStepIdxList_dev[sa16d],
             &reactantSpcIdxList_dev[sa16d],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             destroyOut_dev, scatterAddStreams[3]);
          perf_net_scatterAdd_gpu_atomic_global_8op
            (sa4d-sa8d,
             &reactantStepIdxList_dev[sa8d],
             &reactantSpcIdxList_dev[sa8d],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             destroyOut_dev, scatterAddStreams[4]);
           perf_net_scatterAdd_gpu_atomic_global_4op
            (sa2d-sa4d,
             &reactantStepIdxList_dev[sa4d],
             &reactantSpcIdxList_dev[sa4d],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             destroyOut_dev, scatterAddStreams[5]);
          perf_net_scatterAdd_gpu_atomic_global_2op
            (sa1d-sa2d,
             &reactantStepIdxList_dev[sa2d],
             &reactantSpcIdxList_dev[sa2d],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             destroyOut_dev, scatterAddStreams[6]);
          perf_net_scatterAdd_gpu_atomic_global
            (totReac-sa1d,
             &reactantStepIdxList_dev[sa1d],
             &reactantSpcIdxList_dev[sa1d],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             destroyOut_dev, scatterAddStreams[7]);

          perf_net_scatterAdd_gpu_atomic_global_128op
            (sa64c,
             &productStepIdxList_dev[0],
             &productSpcIdxList_dev[0],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             createOut_dev, scatterAddStreams[8]);
          perf_net_scatterAdd_gpu_atomic_global_64op
            (sa32c-sa64c,
             &productStepIdxList_dev[sa64c],
             &productSpcIdxList_dev[sa64c],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             createOut_dev, scatterAddStreams[9]);
          perf_net_scatterAdd_gpu_atomic_global_32op
            (sa16c-sa32c,
             &productStepIdxList_dev[sa32c],
             &productSpcIdxList_dev[sa32c],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             createOut_dev, scatterAddStreams[10]);
          perf_net_scatterAdd_gpu_atomic_global_16op
            (sa8c-sa16c,
             &productStepIdxList_dev[sa16c],
             &productSpcIdxList_dev[sa16c],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             createOut_dev, scatterAddStreams[11]);
          perf_net_scatterAdd_gpu_atomic_global_8op
            (sa4c-sa8c,
             &productStepIdxList_dev[sa8c],
             &productSpcIdxList_dev[sa8c],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             createOut_dev, scatterAddStreams[12]);
           perf_net_scatterAdd_gpu_atomic_global_4op
            (sa2c-sa4c,
             &productStepIdxList_dev[sa4c],
             &productSpcIdxList_dev[sa4c],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             createOut_dev, scatterAddStreams[13]);
          perf_net_scatterAdd_gpu_atomic_global_2op
            (sa1c-sa2c,
             &productStepIdxList_dev[sa2c],
             &productSpcIdxList_dev[sa2c],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             createOut_dev, scatterAddStreams[14]);
          perf_net_scatterAdd_gpu_atomic_global
            (totProd-sa1c,
             &productStepIdxList_dev[sa1c],
             &productSpcIdxList_dev[sa1c],
             nReactors,
             nStep,
             stepOut_dev,
             nSpc,
             createOut_dev, scatterAddStreams[15]);

    cudaEventRecord(endProd);

    cudaEventRecord(startTxProd);
    if(outputAll)
    {
        cudaMemcpyAsync(createOut_tp,createOut_dev,cpySpcSize,cudaMemcpyDeviceToHost);
        cudaMemcpyAsync(destroyOut_tp,destroyOut_dev,cpySpcSize,cudaMemcpyDeviceToHost);
    }
    cudaEventRecord(endTxProd);

    cudaEventRecord(startNet);
    perf_net_cuda_net_rates(nSpc*nReactors,createOut_dev,destroyOut_dev,netOut_dev);
    cudaEventRecord(endNet);

    cudaEventRecord(startTxNet);
    cudaMemcpy(netOut_tp,netOut_dev,cpySpcSize,cudaMemcpyDeviceToHost);
    cudaEventRecord(endTxNet);

    cudaEventSynchronize(endTxStep);
    cudaEventSynchronize(endTxProd);
    cudaEventSynchronize(endTxNet);
    if(transposeInOut)
    {    
        for(j=0; j < nReactors; ++j)
        {
            for(k=0;k<nSpc;++k)
            {
                int origIdx = (j*nSpc+k);
                int transIdx = (k*nReactors+j);
                netOut[origIdx] = netOut_tp[transIdx];
                if(outputAll)
                {
                    createOut[origIdx] = createOut_tp[transIdx];
                    destroyOut[origIdx] = destroyOut_tp[transIdx];
                }
            }
            if(outputAll)
            {
                for(k=0;k<nStep;++k)
                {
                    int origIdx = (j*nStep+k);
                    int transIdx = (k*nReactors+j);
                    stepOut[origIdx] = stepOut_tp[transIdx];
                }
            }
        }
    }    


    cudaEventSynchronize(endProd);
    cudaEventSynchronize(endNet);

    cudaEventElapsedTime(&cudaTime_ms, startTxStep, endTxStep);
    gpuTxTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startProd, endProd);
    gpuProdTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startTxProd, endTxProd);
    gpuTxTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startNet, endNet);
    gpuNetTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startTxNet, endTxNet);
    gpuTxTime += cudaTime_ms/1000.0;


    if(transposeInOut)
    {
        free(netOut_tp);
        if(outputAll)
        {
            free(createOut_tp);
            free(destroyOut_tp);
            free(stepOut_tp);
        }
    }

    cudaEventDestroy(startTxC);
    cudaEventDestroy(endTxC);
    cudaEventDestroy(startTxStep);
    cudaEventDestroy(endTxStep);
    cudaEventDestroy(startK);
    cudaEventDestroy(endK);
    cudaEventDestroy(startStep);
    cudaEventDestroy(endStep);
    cudaEventDestroy(startTxProd);
    cudaEventDestroy(endTxProd);
    cudaEventDestroy(startProd);
    cudaEventDestroy(endProd);
    cudaEventDestroy(startTxNet);
    cudaEventDestroy(endTxNet);
    cudaEventDestroy(startNet);
    cudaEventDestroy(endNet);

#ifdef ZERORK_FULL_DEBUG
    cudaDeviceSynchronize();
    checkCudaError(cudaGetLastError(),"perf_net_cuda::getRatesFromTC_mr(): end of routine");
#endif
}

void perf_net_cuda::calcRatesFromTC_CUDA_mr_dev(const int nReactors, const double T_out_dev[], const double C_out_dev[],
				   double* netOut_out_dev, double* createOut_out_dev,
				   double* destroyOut_out_dev, double* stepOut_out_dev)
{
    float cudaTime_ms;

    perf_net_cuda_setup_memory(nReactors, nSpc, C_dev, C_out_dev,
                               createOut_out_dev, destroyOut_out_dev,
                               scatterAddStreams[0]);

#ifdef ZERORK_FULL_DEBUG
    cudaDeviceSynchronize();
    checkCudaError(cudaGetLastError(),"perf_net_cuda::calcRatesFromTC_CUDA_mr_dev( ... ) {cudaMemCpy}");
#endif

#ifdef ZERORK_CUDA_EVENTS
    cudaEvent_t startTxStep,endTxStep,startK,endK,startStep,endStep,
                startTxProd,endTxProd,startProd,endProd,startTxNet,endTxNet,startNet,endNet;
    cudaEventCreate(&startTxStep);
    cudaEventCreate(&endTxStep);
    cudaEventCreate(&startK);
    cudaEventCreate(&endK);
    cudaEventCreate(&startStep);
    cudaEventCreate(&endStep);
    cudaEventCreate(&startTxProd);
    cudaEventCreate(&endTxProd);
    cudaEventCreate(&startProd);
    cudaEventCreate(&endProd);
    cudaEventCreate(&startTxNet);
    cudaEventCreate(&endTxNet);
    cudaEventCreate(&startNet);
    cudaEventCreate(&endNet);

    cudaEventRecord(startK);
#endif
    static_cast<rate_const_cuda*>(rateConstPtr)->updateK_CUDA_mr(nReactors,T_out_dev,C_out_dev,stepOut_out_dev);
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endK);
#endif

#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(startStep);
#endif
    //N.B. this happens in default stream so we don't need to sync the streams
    perf_net_cuda_rxn_conc_mult_mr(nReactors, nSpc, nStep, maxReactants,
        reactantSpcIdxListUnwrapped_dev, C_dev, stepOut_out_dev);
#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endStep);
#endif

#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(startTxStep);
    cudaEventRecord(endTxStep);
    cudaDeviceSynchronize();

    cudaEventSynchronize(endK);
    cudaEventSynchronize(endStep);

    cudaEventElapsedTime(&cudaTime_ms, startK, endK);
    gpuKTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startStep, endStep);
    gpuStepTime += cudaTime_ms/1000.0;

    cudaEventRecord(startProd);
#endif

#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endStep);
#endif

#ifndef ZERORK_SCATTER_ADD_FUSED
        perf_net_scatterAdd_gpu_atomic_global_128op
          (sa64d,
           &reactantStepIdxList_dev[0],
           &reactantSpcIdxList_dev[0],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           destroyOut_out_dev, scatterAddStreams[0]);
        perf_net_scatterAdd_gpu_atomic_global_64op
          (sa32d-sa64d,
           &reactantStepIdxList_dev[sa64d],
           &reactantSpcIdxList_dev[sa64d],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           destroyOut_out_dev, scatterAddStreams[1]);
        perf_net_scatterAdd_gpu_atomic_global_32op
          (sa16d-sa32d,
           &reactantStepIdxList_dev[sa32d],
           &reactantSpcIdxList_dev[sa32d],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           destroyOut_out_dev, scatterAddStreams[2]);
        perf_net_scatterAdd_gpu_atomic_global_16op
          (sa8d-sa16d,
           &reactantStepIdxList_dev[sa16d],
           &reactantSpcIdxList_dev[sa16d],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           destroyOut_out_dev, scatterAddStreams[3]);
        perf_net_scatterAdd_gpu_atomic_global_8op
          (sa4d-sa8d,
           &reactantStepIdxList_dev[sa8d],
           &reactantSpcIdxList_dev[sa8d],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           destroyOut_out_dev, scatterAddStreams[4]);
         perf_net_scatterAdd_gpu_atomic_global_4op
          (sa2d-sa4d,
           &reactantStepIdxList_dev[sa4d],
           &reactantSpcIdxList_dev[sa4d],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           destroyOut_out_dev, scatterAddStreams[5]);
        perf_net_scatterAdd_gpu_atomic_global_2op
          (sa1d-sa2d,
           &reactantStepIdxList_dev[sa2d],
           &reactantSpcIdxList_dev[sa2d],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           destroyOut_out_dev, scatterAddStreams[6]);
        perf_net_scatterAdd_gpu_atomic_global
          (totReac-sa1d,
           &reactantStepIdxList_dev[sa1d],
           &reactantSpcIdxList_dev[sa1d],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           destroyOut_out_dev, scatterAddStreams[7]);

        perf_net_scatterAdd_gpu_atomic_global_128op
          (sa64c,
           &productStepIdxList_dev[0],
           &productSpcIdxList_dev[0],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           createOut_out_dev, scatterAddStreams[8]);
        perf_net_scatterAdd_gpu_atomic_global_64op
          (sa32c-sa64c,
           &productStepIdxList_dev[sa64c],
           &productSpcIdxList_dev[sa64c],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           createOut_out_dev, scatterAddStreams[9]);
        perf_net_scatterAdd_gpu_atomic_global_32op
          (sa16c-sa32c,
           &productStepIdxList_dev[sa32c],
           &productSpcIdxList_dev[sa32c],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           createOut_out_dev, scatterAddStreams[10]);
        perf_net_scatterAdd_gpu_atomic_global_16op
          (sa8c-sa16c,
           &productStepIdxList_dev[sa16c],
           &productSpcIdxList_dev[sa16c],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           createOut_out_dev, scatterAddStreams[11]);
        perf_net_scatterAdd_gpu_atomic_global_8op
          (sa4c-sa8c,
           &productStepIdxList_dev[sa8c],
           &productSpcIdxList_dev[sa8c],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           createOut_out_dev, scatterAddStreams[12]);
         perf_net_scatterAdd_gpu_atomic_global_4op
          (sa2c-sa4c,
           &productStepIdxList_dev[sa4c],
           &productSpcIdxList_dev[sa4c],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           createOut_out_dev, scatterAddStreams[13]);
        perf_net_scatterAdd_gpu_atomic_global_2op
          (sa1c-sa2c,
           &productStepIdxList_dev[sa2c],
           &productSpcIdxList_dev[sa2c],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           createOut_out_dev, scatterAddStreams[14]);
        perf_net_scatterAdd_gpu_atomic_global
          (totProd-sa1c,
           &productStepIdxList_dev[sa1c],
           &productSpcIdxList_dev[sa1c],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           createOut_out_dev, scatterAddStreams[15]);
#else
        perf_net_scatterAdd_gpu_atomic_global_fused
          (maxOpsDestroy, nOpsDestroy_dev,
           &reactantStepIdxList_dev[0],
           &reactantSpcIdxList_dev[0],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           destroyOut_out_dev, scatterAddStreams[0]);

        perf_net_scatterAdd_gpu_atomic_global_fused
          (maxOpsCreate, nOpsCreate_dev,
           &productStepIdxList_dev[0],
           &productSpcIdxList_dev[0],
           nReactors,
           nStep,
           stepOut_out_dev,
           nSpc,
           createOut_out_dev, scatterAddStreams[8]);
#endif

#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endProd);

    cudaEventRecord(startTxProd);
    cudaEventRecord(endTxProd);

    cudaEventRecord(startNet);
#endif
    //N.B. this happens in default stream so we don't need to sync the streams
    perf_net_cuda_net_rates(nSpc*nReactors, createOut_out_dev,
        destroyOut_out_dev, netOut_out_dev);
#ifdef ZERORK_FULL_DEBUG
    cudaDeviceSynchronize();
    checkCudaError(cudaGetLastError(),"cuda_net_rates (multi_dev)");
#endif

#ifdef ZERORK_CUDA_EVENTS
    cudaEventRecord(endNet);

    cudaEventRecord(startTxNet);
    cudaEventRecord(endTxNet);

    cudaEventSynchronize(endTxStep);
    cudaEventSynchronize(endTxProd);
    cudaEventSynchronize(endTxNet);

    cudaEventSynchronize(endProd);
    cudaEventSynchronize(endNet);

    cudaEventElapsedTime(&cudaTime_ms, startTxStep, endTxStep);
    gpuTxTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startProd, endProd);
    gpuProdTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startTxProd, endTxProd);
    gpuTxTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startNet, endNet);
    gpuNetTime += cudaTime_ms/1000.0;
    cudaEventElapsedTime(&cudaTime_ms, startTxNet, endTxNet);
    gpuTxTime += cudaTime_ms/1000.0;

    cudaEventDestroy(startTxStep);
    cudaEventDestroy(endTxStep);
    cudaEventDestroy(startK);
    cudaEventDestroy(endK);
    cudaEventDestroy(startStep);
    cudaEventDestroy(endStep);
    cudaEventDestroy(startTxProd);
    cudaEventDestroy(endTxProd);
    cudaEventDestroy(startProd);
    cudaEventDestroy(endProd);
    cudaEventDestroy(startTxNet);
    cudaEventDestroy(endTxNet);
    cudaEventDestroy(startNet);
    cudaEventDestroy(endNet);
#endif

    //N.B. : This seems like it is pointless, but in testing
    //it improves performance.  Hypothesis is that it keeps
    //the cpu from getting to far ahead
    cudaStreamSynchronize((cudaStream_t) 0);
#ifdef ZERORK_FULL_DEBUG
    cudaDeviceSynchronize();
    checkCudaError(cudaGetLastError(),"perf_net_cuda::getRatesFromTC_mr_dev(): end of routine");
#endif
}


typedef struct
{
  int x_index;
  int y_index;
} indexPair;


// int compareXthenY(void *a, void *b);
//
// A comparison function supplied to qsort() that compares the index pairs
// 'a' and 'b' and returns -1, 0 or 1 if 'a' is less than, equal to, or 
// greater than 'b'.  The comparison is first done on the x index and then
// the y index.
int compareXthenY(const void *a, const void *b)
{
  indexPair *aptr = (indexPair *)a;
  indexPair *bptr = (indexPair *)b;
  if(aptr->x_index < bptr->x_index)
    {return -1;}  // a < b
  if(aptr->x_index > bptr->x_index)
    {return 1;}   // a < b

  // x_index is equal, now compare y_index
  if(aptr->y_index < bptr->y_index)
    {return -1;}  // a < b
  if(aptr->y_index > bptr->y_index)
    {return 1;}   // a < b

  return 0;
}

// int compareYthenX(void *a, void *b);
//
// A comparison function supplied to qsort() that compares the index pairs
// 'a' and 'b' and returns -1, 0 or 1 if 'a' is less than, equal to, or 
// greater than 'b'.  The comparison is first done on the y index and then
// the x index.
static int compareYthenX(const void *a, const void *b)
{
  indexPair *aptr = (indexPair *)a;
  indexPair *bptr = (indexPair *)b;
  if(aptr->y_index < bptr->y_index)
    {return -1;}  // a < b
  if(aptr->y_index > bptr->y_index)
    {return 1;}   // a < b

  // y_index is equal, now compare x_index
  if(aptr->x_index < bptr->x_index)
    {return -1;}  // a < b
  if(aptr->x_index > bptr->x_index)
    {return 1;}   // a < b

  return 0;
}
// void indexPairSort(const int len,
//                    const char firstComp,
//                    int x[],
//                    int y[])
//
// Sort the index lists x and y as pairs sorting them first by x or y
// depending on the character string provided as firstComp.  If firstComp is
// not equal to 'y' or 'Y', then the x-first sort will be used.
static void indexPairSort(const int len,
		   const char firstComp,
		   int x[],
		   int y[])
{
  int j;
  std::vector<indexPair> tmpList(len);

  // store the values of x and y into the index pair list
  for(j=0; j<len; j++)
    {
      tmpList[j].x_index=x[j];
      tmpList[j].y_index=y[j];      
    }

  // perform a quick sort on the index pairs
  if(firstComp == 'y' || firstComp == 'Y') 
    { // sort by the y_index first
      qsort((void *)&tmpList[0],
	    len,
	    sizeof(indexPair),
	    compareYthenX);
    }
  else
    { // sort by the x_index first (default)
      qsort((void *)&tmpList[0],
	    len,
	    sizeof(indexPair),
	    compareXthenY);
    }

  // copy the sorted result back in the original x and y arrays
  for(j=0; j<len; j++)
    {
      x[j]=tmpList[j].x_index;
      y[j]=tmpList[j].y_index;      
    }

}

// re-order the source and destination index pairs in n-tuples by identical
// destinations.  Destination addresses with operation counts with remainders
// are aggregated at the end. The starting address in the operation list is
// returned.
int perf_net_cuda::reorderScatterAdd_by_ntuple(const int ntuple,
			     const int nOps,
			     const int srcSize,
			     const int destSize,
			     int srcId[],
			     int destId[])
{
  int j,k;
  int remId,groupId,nRem,nGroup;
  int ntupleCount;
  int *destCount;
  int *newSrcId,*newDestId;

  destCount = (int *)malloc(sizeof(int)*destSize);
  newSrcId  = (int *)malloc(sizeof(int)*nOps);
  newDestId = (int *)malloc(sizeof(int)*nOps);

  // sort the oepration pairs (destId, srcId) by the destId first
  indexPairSort(nOps,'x',&destId[0],&srcId[0]);

  // initialize the counter array for the destination indexes
  for(j=0; j<destSize; j++)
    {destCount[j]=0;}
  
  // count the distribution of destination indexes
  for(j=0; j<nOps; j++)
    {destCount[destId[j]]++;}

  // count the number of complete ntuples
  ntupleCount=0;
  for(j=0; j<destSize; j++)
    {ntupleCount+=(destCount[j]/ntuple);}
  
  // regroup the sorted index pairs so that all the complete n-tuples
  // are located at the beginning and all the partial n-tuple remainders
  // are at the end
  j = 0;                      // number of operation index pairs places
  groupId = 0;                // starting address of the n-tuple groups
  remId = ntuple*ntupleCount; // starting address of the remainders
  while(j < nOps)
    {
      // nGroup is the number of index pairs in full n-tuples
      // nRem   is the number of leftover index pairs that fail to fill
      //        an n-tuple
      nGroup  = destCount[destId[j]]/ntuple;
      nGroup *= ntuple;
      nRem    = destCount[destId[j]]-nGroup;

      // place the n-tuple groups
      for(k=0; k<nGroup; k++) 
	{
	  newDestId[groupId] = destId[j];
	  newSrcId[groupId]  = srcId[j];
	  ++j; ++groupId;
	}
      // place the remainders
      for(k=0; k<nRem; k++)
	{
	  newDestId[remId] = destId[j];
	  newSrcId[remId]  = srcId[j];
	  ++j; ++remId;
	}
    }

  // copy the regrouped index pairsback to the 
  for(j=0; j<nOps; j++)
    {
      destId[j] = newDestId[j];
      srcId[j]  = newSrcId[j];
    }

  //  for(j=0; j<nOps; j++)
  //{
  //  printf("op list %5d: (dest = %5d, src = %5d)\n",
  // 	     j,destId[j],srcId[j]);
  //} 
//  printf("# %d-tuple fraction: %8.6f (count %d)\n",
//	 ntuple,
//	 (double)(ntuple*ntupleCount)/(double)nOps,
//	 ntuple*ntupleCount);
	 
    
  free(newDestId);
  free(newSrcId);
  free(destCount);

  return ntuple*ntupleCount;
}


} // namespace zerork

