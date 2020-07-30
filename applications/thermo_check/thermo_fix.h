#ifndef THERMO_FIX_H_
#define THERMO_FIX_H_

double SpecificHeatDelta(const double T,
                         const double coef_low[],
                         const double coef_high[]);
double SpecificHeatDerivativeDelta(const double T,
                                   const double coef_low[],
                                   const double coef_high[]);
double EntropyDelta(const double T,
                    const double coef_low[],
                    const double coef_high[]);

void RefitKeepHigh(const int num_points,
                   const int T_resolution,
                   const double T_fixed,
                   const double T_min,
                   const double T_match,
                   const double T_max,
                   const double coef_low[],
                   const double coef_high[],
                   double *refit_T_match,
                   double refit_coef_low[],
                   double refit_coef_high[]);
void RefitKeepHigh(const int num_points,
                   const int T_resolution,
                   const double T_fixed,
                   const JanafThermoData *original,
                   const bool find_new_match,
		   JanafThermoData *refit);
void RefitKeepHighGlobalTMatch(const int num_points,
                               const int T_resolution,
                               const double T_fixed,
                               const JanafThermoData *original,
                               const double global_T_match,
		               JanafThermoData *refit);

void RefitKeepHighGlobalTMatch(const int num_points,
                               const int T_resolution,
                               const double T_fixed,
                               const double T_min,
                               const double T_match,
                               const double T_max,
                               const double coef_low[],
                               const double coef_high[],
                               const double global_T_match,
                               double refit_coef_low[],
                               double refit_coef_high[]);

void MatchEnthalpy(const double T_fixed,
                   const double H_fixed,
                   double coef[]);
void MatchEntropy(const double T_fixed,
                  const double S_fixed,
                  double coef[]);
double MinSpecificHeatDelta(const double T_tol,
                            const double T_min,
                            const double T_match,
                            const double T_max,
                            const double coef_low[],
                            const double coef_high[]);
bool PolynomialFit(const int num_points,
                   const int degree,
                   const int num_fixed_points,
                   const int num_fixed_slopes,
                   const double x[],
                   const double f[],
                   const double x_fixed_points[],
                   const double f_fixed_points[],
                   const double x_fixed_slopes[],
                   const double dfdx_fixed_slopes[],
                   double coef[]);

double RoundToTenPower(const int n, const double a);


#endif
