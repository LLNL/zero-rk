#ifdef _WIN32
#define _USE_MATH_DEFINES //for M_PI
#endif
#include <math.h>

#include "zerork/constants.h"

#include "collision_integrals.h"

namespace transport {

double getOmega2_2_LJ(const double at_reduced_temperature)
{
  int range_id;
  double range_fraction,log_omega,log_omega_min;
  // find the reduced temperature range for interpolation (or extrapolation)
  // range_id is set so that the reduced_temperature is found between points
  // range_id and range_id+1.  This search is performed as a brute force
  // calculation.
  //
  // TODO - use a more efficient search/data structure.
  if(at_reduced_temperature < KReducedTemperature[1]) {
    range_id = 0; 
  } else if(at_reduced_temperature >= KReducedTemperature[KNumOmegaPoints-2]) {
    range_id = KNumOmegaPoints-2;
  } else {
    for(int j=1; j<KNumOmegaPoints-2; ++j) {
      if(KReducedTemperature[j] <= at_reduced_temperature &&
         at_reduced_temperature < KReducedTemperature[j+1]) {
        range_id = j;
        break;
      } 
    }
  } // end of else (brute force search)

  // TODO - storing the logarithms of the reduced temperature and Omega_2_2
  // constants will further speed up computation.
  log_omega_min  = log(KOmega2_2_LJ[range_id]);
  range_fraction = log(at_reduced_temperature/KReducedTemperature[range_id])/
    log(KReducedTemperature[range_id+1]/KReducedTemperature[range_id]);

  log_omega = log_omega_min + range_fraction*log(KOmega2_2_LJ[range_id+1]/
                                                 KOmega2_2_LJ[range_id]);
  return exp(log_omega);
}

// Returns the viscosity for the Lennard-Jones 6-12 potential using the
// tabulated Omega(2,2) collision integrals.
//
// Inputs:  temperature [K]
//          molecular weight [kg/kmol]
//          epsilon/k_boltzmann [K]
//          sigma [m] (note that sigma is often reported in Angstroms)
// Outputs: viscosity [Pa*s]
double getViscosity_LJ(const double temperature,
                       const double molecular_weight,
                       const double epsilon_k,
                       const double sigma)
{
  const double conversion_factor=
    5.0/16.0*sqrt(zerork::KBoltzmann/zerork::KAvogadroNumber/M_PI);
  double collision_term = getOmega2_2_LJ(temperature/epsilon_k);

  return conversion_factor*sqrt(molecular_weight*temperature)/
    (collision_term*sigma*sigma);
}

} // end of namespace transport
