#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "polynomial.h"
#include "janaf_thermo.h"

// Evaluate the non-dimensional specific heat at constant pressure from the 
// 4th degree polynomial defined by coef[].
//
// Cp/R = coef[0] + coef[1]*T + coef[2]*T**2 + coef[3]*T**3 + coef[4]*T**4
double SpecificHeat(const double T, const double coef[])
{
  return coef[0] + T*(coef[1] + T*(coef[2] + T*(coef[3] + T*coef[4])));
}
// Evaluate the temperature derivative of the non-dimensional specific heat 
// at constant pressure from the 4th degree polynomial defined by coef[].
//
// Cp/R = coef[0] + coef[1]*T + coef[2]*T**2 + coef[3]*T**3 + coef[4]*T**4
double SpecificHeatDerivative(const double T, const double coef[])
{
  return coef[1] + T*(coef[2]*2.0 + T*(coef[3]*3.0 + T*coef[4]*4.0));
}
//Evaluate the enthalpy H/RT (normalized by the gas constant and temperature)
//at a given temperature T and single set of JANAF coefficients
double Enthalpy(const double T, const double coef[])
{
  const double ONE_THIRD = 1.0/3.0;
  
  return coef[0] + coef[5]/T + T*(coef[1]*0.5 + 
                                  T*(coef[2]*ONE_THIRD +
                                     T*(coef[3]*0.25 + 
                                        T*coef[4]*0.2)));
}
// Evaluate the entropy S/R (normalized by the gas constant) at a given 
// temperature T and single set of JANAF coefficients
double Entropy(const double T, const double coef[])
{
  const double ONE_THIRD = 1.0/3.0;
  
  return coef[0]*log(T) + coef[6] + T*(coef[1] + 
                                       T*(coef[2]*0.5 +
                                          T*(coef[3]*ONE_THIRD +
                                             T*coef[4]*0.25)));
}

double JanafSpecificHeat(const double T, const JanafThermoData *data)
{
  if(T < data->T_match) {
    return SpecificHeat(T,data->low_coef);
  }
  return SpecificHeat(T,data->high_coef);
}
double JanafEnthalpy(const double T, const JanafThermoData *data)
{
  if(T < data->T_match) {
    return Enthalpy(T,data->low_coef);
  }
  return Enthalpy(T,data->high_coef);
}
double JanafEntropy(const double T, const JanafThermoData *data)
{
  if(T < data->T_match) {
    return Entropy(T,data->low_coef);
  }
  return Entropy(T,data->high_coef);
}

double JanafSpecificHeatClipped(const double T, const JanafThermoData *data)
{
  if(T < data->T_min || T > data->T_max) return NAN;
  if(T < data->T_match) {
    return SpecificHeat(T,data->low_coef);
  }
  return SpecificHeat(T,data->high_coef);
}
double JanafEnthalpyClipped(const double T, const JanafThermoData *data)
{
  if(T < data->T_min || T > data->T_max) return NAN;
  if(T < data->T_match) {
    return Enthalpy(T,data->low_coef);
  }
  return Enthalpy(T,data->high_coef);
}
double JanafEntropyClipped(const double T, const JanafThermoData *data)
{
  if(T < data->T_min || T > data->T_max) return NAN;
  if(T < data->T_match) {
    return Entropy(T,data->low_coef);
  }
  return Entropy(T,data->high_coef);
}

double MaxDifference(double (*janaf_func)(const double, 
                                          const JanafThermoData *),
                     const int num_samples,
                     const JanafThermoData *data1,
                     const JanafThermoData *data2,
                     double *T_max_difference)
{
  double func_max_difference, T_current, func1_eval,func2_eval,func_difference;
  double T_min,T_max;

  func_max_difference = -1.0e300;

  T_min = ((data1->T_min < data2->T_min) ? data1->T_min : data2->T_min);
  T_max = ((data1->T_max > data2->T_max) ? data1->T_max : data2->T_max);

  for(int j=0; j<num_samples; ++j) {

    T_current = T_min + (T_max-T_min)*static_cast<double>(j)/(num_samples-1.0);
    func1_eval = (*janaf_func)(T_current,data1);
    func2_eval = (*janaf_func)(T_current,data2);

    func_difference = fabs(func1_eval-func2_eval);
    if(func_difference > func_max_difference) {
      func_max_difference = func_difference;
      *T_max_difference = T_current;
    }
  }
  return func_max_difference;
}

void PrintAllThermoFunctions(const int num_points,
                             const double T_min,
                             const double T_match,
                             const double T_max,
                             const double coef_low[],
                             const double coef_high[])
{
  double T;
  double Cp_R, H_RT, S_R, dCp_dT_R;
  printf("# lower temperature range\n");
  for(int j=0; j<num_points; ++j) {
    T = T_min + (T_match-T_min)*static_cast<double>(j)/(num_points-1.0);
    Cp_R     = SpecificHeat(T,coef_low);
    H_RT     = Enthalpy(T,coef_low);
    S_R      = Entropy(T,coef_low);
    dCp_dT_R = SpecificHeatDerivative(T,coef_low);
    printf("%12.6f  %12.8f  %12.8f  %12.8f  %16.8e\n",
           T,
           Cp_R,
           H_RT,
           S_R,
           dCp_dT_R);
  }
  printf("# higher temperature range\n");
  for(int j=0; j<num_points; ++j) {
    T = T_match + (T_max-T_match)*static_cast<double>(j)/(num_points-1.0);
    Cp_R     = SpecificHeat(T,coef_high);
    H_RT     = Enthalpy(T,coef_high);
    S_R      = Entropy(T,coef_high);
    dCp_dT_R = SpecificHeatDerivative(T,coef_high);
    printf("%12.6f  %12.8f  %12.8f  %12.8f  %16.8e\n",
           T,
           Cp_R,
           H_RT,
           S_R,
           dCp_dT_R);
  }
}

void PrintJanaf(const JanafThermoData &data)
{
  printf("JANAF data:\n");
  printf("T_min   = %.18g\n",data.T_min);
  printf("T_match = %.18g\n",data.T_match);
  printf("T_max   = %.18g\n",data.T_max);
  printf("low_coef[0]  = %.18g\n",data.low_coef[0]);
  printf("low_coef[1]  = %.18g\n",data.low_coef[1]);
  printf("low_coef[2]  = %.18g\n",data.low_coef[2]);
  printf("low_coef[3]  = %.18g\n",data.low_coef[3]);
  printf("low_coef[4]  = %.18g\n",data.low_coef[4]);
  printf("low_coef[5]  = %.18g\n",data.low_coef[5]);
  printf("low_coef[6]  = %.18g\n",data.low_coef[6]);
  printf("high_coef[0]  = %.18g\n",data.high_coef[0]);
  printf("high_coef[1]  = %.18g\n",data.high_coef[1]);
  printf("high_coef[2]  = %.18g\n",data.high_coef[2]);
  printf("high_coef[3]  = %.18g\n",data.high_coef[3]);
  printf("high_coef[4]  = %.18g\n",data.high_coef[4]);
  printf("high_coef[5]  = %.18g\n",data.high_coef[5]);
  printf("high_coef[6]  = %.18g\n",data.high_coef[6]);
  fflush(stdout);
}

// Compute the jump from low to high (high - low) at T_match for
// Cp/R = jump[0]
// H/RT = jump[1]
// S/R  = jump[2]
void JanafJump(const JanafThermoData &data, double jump[])
{
  jump[0] = SpecificHeat(data.T_match,data.high_coef) -
    SpecificHeat(data.T_match,data.low_coef);
  jump[1] = Enthalpy(data.T_match,data.high_coef) -
    Enthalpy(data.T_match,data.low_coef);
  jump[2] = Entropy(data.T_match,data.high_coef) -
    Entropy(data.T_match,data.low_coef);
                         

}
// Compute the number of extrema found in the Cp/R polynomial fit
// for a give JANAF definition. Store the temperature and specific heats
// for any extrema found.
int GetJanafExtrema(const JanafThermoData &data,
                    const double imag_root_rtol,
                    const double imag_root_atol,
                    std::vector<double> *temperatures,
                    std::vector<double> *specific_heats)
{
  const int max_degree = 4;
  double min_specific_heat,max_specific_heat;
  int num_extrema;
  double temperature_extrema[max_degree-1];
  double specific_heat_extrema[max_degree-1];

  // clear vectors storing the extrema
  temperatures->clear();
  specific_heats->clear();

  // -------------------------------------------------------------------------
  // low temperature branch
  min_specific_heat = SpecificHeat(data.T_min,data.low_coef);
  if(min_specific_heat <= 0.0) {
    // record negative specific heat at T_min as an extrema
    temperatures->push_back(data.T_min);
    specific_heats->push_back(min_specific_heat);
  }
  num_extrema = GetPolynomialExtrema(max_degree,
                                     data.low_coef,
                                     imag_root_rtol,
                                     imag_root_atol,
                                     temperature_extrema,
                                     specific_heat_extrema);
  if(num_extrema < 0) {
    printf("WARNING: In GetJanafExtrema(...),\n");
    printf("         GetPolynomialExtrema(...) returned %d\n", num_extrema);
    printf("         for number of extrema for the low temperature range.\n");
    printf("         Ignoring interior extrema in this region.\n");
    fflush(stdout);
  }
  for(int j=0; j<num_extrema; ++j) {
    if(data.T_min < temperature_extrema[j] &&
       temperature_extrema[j] < data.T_match) {
      // record interior extrema
      temperatures->push_back(temperature_extrema[j]);
      specific_heats->push_back(specific_heat_extrema[j]);
    }
  }
  // max refers to the maximum of the temperature range
  max_specific_heat = SpecificHeat(data.T_match,data.low_coef);
  if(max_specific_heat < min_specific_heat) {
    // record monotonically decreasing end point as an extrema
    temperatures->push_back(data.T_match);
    specific_heats->push_back(max_specific_heat);  
  }
  // -------------------------------------------------------------------------
  // high temperature branch
  min_specific_heat = SpecificHeat(data.T_match,data.high_coef);
  if(min_specific_heat <= 0.0) {
    // record negative specific heat at T_min as an extrema
    temperatures->push_back(data.T_match);
    specific_heats->push_back(min_specific_heat);
  }
  num_extrema = GetPolynomialExtrema(max_degree,
                                     data.high_coef,
                                     imag_root_rtol,
                                     imag_root_atol,
                                     temperature_extrema,
                                     specific_heat_extrema);
  if(num_extrema < 0) {
    printf("WARNING: In GetJanafExtrema(...),\n");
    printf("         GetPolynomialExtrema(...) returned %d\n", num_extrema);
    printf("         for number of extrema for the low temperature range.\n");
    printf("         Ignoring interior extrema in this region.\n");
    fflush(stdout);
  }
  for(int j=0; j<num_extrema; ++j) {
    if(data.T_match < temperature_extrema[j] &&
       temperature_extrema[j] < data.T_max) {
      // record interior extrema
      temperatures->push_back(temperature_extrema[j]);
      specific_heats->push_back(specific_heat_extrema[j]);
    }
  }
  // max refers to the maximum of the temperature range
  max_specific_heat = SpecificHeat(data.T_max,data.high_coef);
  if(max_specific_heat < min_specific_heat) {
    // record monotonically decreasing end point as an extrema
    temperatures->push_back(data.T_max);
    specific_heats->push_back(max_specific_heat);  
  }

  return static_cast<int>(temperatures->size());
}
