#ifndef ZERORK_MECHANISM_KERNELS_H
#define ZERORK_MECHANISM_KERNELS_H

namespace zerork {

void minusEqualOne(const int n, double *A_dev);
void invert(const int n, double *A_dev);
void meanCpMass_mr(const int nReactors, const int nSpc, const double *RuInvMolWt_dev,
    const double *y_dev, double *cpSpc_dev, double *cpReactors_dev);
void Cv_from_Cp_mr(const int nReactors, const int nSpc, const double *RuInvMolWt_dev,
    const double *y_dev, double *cvSpc_dev, double *cvReactors_dev);
void getCfromVY_mr_dev_wrapper(const int nReactors, const int nSpc,
                       const double *v_dev, const double *y_dev,
                       const double *invMolWt_dev, double *c_dev);
void getMassIntEnergyFromTY_mr_dev_wrapper(const int nReactors, const int nSpc,
                       const double *T_dev, const double *y_dev,
                       const double *RuInvMolWt_dev, double* u_spc_dev, double *u_dev);
void getMassEnthalpyFromTY_mr_dev_wrapper(const int nReactors, const int nSpc,
                       const double *T_dev, const double *y_dev,
                       const double *RuInvMolWt_dev, double* h_spc_dev, double *h_dev);
void getMolWtMixFromY_mr_dev_wrapper(const int nReactors, const int nSpc,
                       const double *y_dev, const double *invMolWt_dev,
                       double *mwMix_dev);
void getDensityfromTPMW_wrapper(const int nReactors, const double Ru,
                                const double T_dev[],
                                const double P_dev[], double dens_dev[]);

void getTemperatureFromEY_mr_dev_part1_wrapper(const int nReactors,
                                               const double T_min, const double T_max,
                                               const double* E,
                                               const double* Emin, const double* Emax,
                                               const double* cv_min, const double* cv_max,
                                               double* T, int* converged);

void getTemperatureFromEY_mr_dev_iter_wrapper(const int nReactors, double tolerance,
                                              const double* E, const double* Eiter,
                                              const double* cv_iter,
                                              double* T, int* converged);

void getPressureFromTVY_wrapper(const int nReactors, const double Ru,
                                const double* T, const double* v, double* P);

} // namespace zerork
#endif
