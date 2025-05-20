#ifndef ZERORK_MECHANISM_CUDA_H
#define ZERORK_MECHANISM_CUDA_H

#include "mechanism.h"

namespace zerork {

class mechanism_cuda : public mechanism
{
 public:
  mechanism_cuda(const char *mechFileName,
		const char *thermFileName, 
		const char *convertFileName,
                int verbosity,
                int nReactorsMax);
  virtual ~mechanism_cuda();

  //Multi-reactor code.
  void getReactionRates_CUDA_mr_dev(const int nReactors, const double T_dev[], const double C_dev[],
             double netOut_dev[], double createOut_dev[],
             double destroyOut_dev[], double stepOut_dev[]);
  void getReactionRatesLimiter_CUDA_mr_dev(const int nReactors, const double T_dev[], const double C_dev[], const double stepLimiter_dev[],
             double netOut_dev[], double createOut_dev[],
             double destroyOut_dev[], double stepOut_dev[]);
  void getMassCpFromTY_mr_dev(const int nReactors, const double T_dev[], const double y_dev[],
       double cvSpc_dev[], double cvReactors_dev[]) const;
  void getMassCvFromTY_mr_dev(const int nReactors, const double T_dev[], const double y_dev[],
       double cvSpc_dev[], double cvReactors_dev[]) const;
  void getEnthalpy_RT_mr_dev(const int nReactors, const double T_dev[], double h_RT_dev[]) const;
  void getIntEnergy_RT_mr_dev(const int nReactors, const double T_dev[], double u_RT_dev[]) const;
  void getCfromVY_mr_dev(const int nReactors, const double v_dev[], const double y_dev[], double c_dev[]) const;
  void getMassIntEnergyFromTY_mr_dev(const int nReactors, const double T_dev[], const double y_dev[],
                                     double u_spc_dev[], double u_dev[] ) const;
  void getMassEnthalpyFromTY_mr_dev(const int nReactors, const double T_dev[], const double y_dev[],
                                    double h_spc_dev[], double h_dev[] ) const;
  void getMolWtMixFromY_mr_dev(const int nReactors, const double y_dev[], double molWtMix_dev[]) const;
  void getDensityFromTPY_mr_dev(const int nReactors, const double T_dev[],
                                const double P_dev[],const double y_dev[],
                                double dens_dev[]) const;
  void getTemperatureFromEY_mr_dev(const int nReactors, const double E[],
                                   const double y[], double temperature_out[]) const;
  void getTemperatureFromHY_mr_dev(const int nReactors, const double H[],
                                   const double y[], double temperature_out[]) const;
  void getPressureFromTVY_mr_dev(const int nReactors, const double* T, const double* v,
                                               const double* y, double* P) const;
  int nReactorsMax();

 private:
  void allocate_device_mem();
  void free_device_mem();

  // CUDA device arrays
  double *molWt_dev;
  double *invMolWt_dev;
  double *RuInvMolWt_dev;

  virtual void initialize_ptrs_cuda(ckr::CKReader *ckrobj, int nReactorsMax);
};

} // namespace zerork

#endif
