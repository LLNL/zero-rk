#ifndef ZERORK_PERF_NET_H
#define ZERORK_PERF_NET_H

#include "info_net.h"
#include "rate_const.h"
#include "external_funcs.h"

namespace zerork {


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

  // Compute the reaction rates and the rate-of-progress of each step after
  // applying the step limiter to rate coefficient.  The step limiter is
  // applied as follows (K is the rate coefficient, K_lim is the limited rate
  // coefficient used to compute the rate of progress of each step, and L is
  // the limit value stored in the step_limiter array):
  //
  //    K_lim = K * [L/(K + L)]  note K_lim <= L
  //
  // the relative error using the limiter is (K/L)/((K/L) + 1) so when
  // K is 1/100 of L the relative error is 0.01/1.01 or ~1%
  //
  void calcRatesFromTC_StepLimiter(const double T,
                                   const double C[],
		                   const double step_limiter[],
			           double netOut[],
                                   double createOut[],
			           double destroyOut[],
                                   double stepOut[]);

  void calcRatesFromTC_StepLimiter_perturbROP(const double T,
                                              const double C[],
                                              const double step_limiter[],
                                              const double perturbMult[],
                                              double netOut[],
                                              double createOut[],
                                              double destroyOut[],
                                              double stepOut[]);


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
  int niTotProd;
  int niTotReac;

  // index storage arrays
  int *reactantSpcIdxList;
  int *reactantStepIdxList;
  int *productSpcIdxList;
  int *productStepIdxList;

  std::vector<int> niReactantSpcIdxList;
  std::vector<int> niReactantStepIdxList;
  std::vector<double> niReactantStoichNumList;
  std::vector<int> niProductSpcIdxList;
  std::vector<int> niProductStepIdxList;
  std::vector<double> niProductStoichNumList;

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

  // special reaction handling
  bool use_non_integer_network_;
  NonIntegerReactionNetwork non_integer_network_;
};

} // namespace zerork

#endif
