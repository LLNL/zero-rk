#ifndef JANAF_THERMO_H_
#define JANAF_THERMO_H_

#include <vector>

typedef struct
{
  double T_min;
  double T_match;
  double T_max;
  double low_coef[7];
  double high_coef[7];
} JanafThermoData;


double SpecificHeat(const double T, const double coef[]);
double SpecificHeatDerivative(const double T, const double coef[]);
double Enthalpy(const double T, const double coef[]);
double Entropy(const double T, const double coef[]);

double JanafSpecificHeat(const double T, const JanafThermoData *data);
double JanafEnthalpy(const double T, const JanafThermoData *data);
double JanafEntropy(const double T, const JanafThermoData *data);

double JanafSpecificHeatClipped(const double T, const JanafThermoData *data);
double JanafEnthalpyClipped(const double T, const JanafThermoData *data);
double JanafEntropyClipped(const double T, const JanafThermoData *data);

void PrintJanaf(const JanafThermoData &data);

double MaxDifference(double (*janaf_func)(const double, 
                                          const JanafThermoData *),
                     const int num_samples,
                     const JanafThermoData *data1,
                     const JanafThermoData *data2,
                     double *T_max_difference);
void JanafJump(const JanafThermoData &data, double jump[]);

int GetJanafExtrema(const JanafThermoData &data,
                    const double imag_root_rtol,
                    const double imag_root_atol,
                    std::vector<double> *temperatures,
                    std::vector<double> *specific_heats);

void PrintAllThermoFunctions(const int num_points,
                             const double T_min,
                             const double T_match,
                             const double T_max,
                             const double coef_low[],
                             const double coef_high[]);



#endif
