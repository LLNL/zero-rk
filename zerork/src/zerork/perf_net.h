#ifndef ZERORK_PERF_NET_H
#define ZERORK_PERF_NET_H

#include "info_net.h"
#include "rate_const.h"
#include "external_funcs.h"

namespace zerork {

typedef struct
{
  int reactantIdx;
  int stepIdx;
  int specIdx;
} multSort;

int compare_multSort(const void *x, const void *y);
int compareSpecStep_multSort(const void *x, const void *y);
int compareStepSpec_multSort(const void *x, const void *y);


class perf_net
{
 public: 
  perf_net(info_net &netobj, rate_const &Kobj);
  virtual ~perf_net();

  void calcRatesFromTC(const double T, const double C[], double netOut[],
		       double createOut[], double destroyOut[],
		       double stepOut[]);
  void calcRatesFromTCM(const double T,
                        const double C[],
                        const double C_mix,
                        double netOut[],
		        double createOut[],
                        double destroyOut[],
		        double stepOut[]);
  void calcRatesFromExplicit(const double T, const double C[], double netOut[],
		       double createOut[], double destroyOut[],
		       double stepOut[]);
  void calcRatesFromTC_unroll16(const double T, const double C[],
           double netOut[], double createOut[], double destroyOut[],
	   double stepOut[]);
  void calcRatesFromTC_perturbROP(const double T, const double C[],
           const double perturbMult[], double netOut[], double createOut[],
           double destroyOut[], double stepOut[]);

//  void writeExplicitRateFunc(const char *fileName, const char *funcName);
//  void writeExplicitRateFunc_minAssign(const char *fileName,
//				       const char *funcName);
  void writeExternalFuncs();
  void write_func_check(FILE* fptr);
  void write_func_rates(FILE* fptr);

  void setUseExRates(){ use_external_rates = true; };
  void unsetUseExRates(){ use_external_rates = false; };
  void setExRatesFunc(external_func_rates_t fn_handle) { ex_func_calc_rates = fn_handle; };

 protected:
  int nStep;
  int nSpc;
  int totProd;
  int totReac;
  int maxReactants;
  
  // index storage arrays
  int *reactantSpcIdxList;
  int *reactantStepIdxList;
  int *productSpcIdxList;
  int *productStepIdxList;
   

  // working arrays:
  double *stepRate;
  double *netSpcProdRate;
  double *spcDestroyRate;
  double *spcCreateRate;
  
  rate_const *rateConstPtr;
  info_net *infoPtr;
    
  bool use_external_rates;
  external_func_rates_t ex_func_calc_rates;

  double cpuKTime,cpuStepTime,cpuProdTime,cpuNetTime;
};

} // namespace zerork

#endif
