#ifndef COLLISION_INTEGRALS_H_
#define COLLISION_INTEGRALS_H_

namespace transport {

double getOmega2_2_LJ(const double reduced_temperature);
double getViscosity_LJ(const double temperature,
                       const double molecular_weight,
                       const double epsilon_k,
                       const double sigma);

// the Omega(2,2) collision integral is the non-dimensional term collision term
// used to compute the first order approximation for viscosity and conductivity
// from Chapman-Enskog theory.  The values below are for a Lennard-Jones 6-12
// potential and correspond to the reduced temperature points in the following
// array.
//
// viscosity = (5/16)*sqrt(pi*m*k*T)/(pi*sigma*sigma*Omega(2,2))
//
//  in [g/(cm-s)] = 2.6693e-5 *sqrt(M*T)/(sigma*sigma*Omega(2,2)),
//                  with M [g/mol], T [K], and sigma [Angstroms]
//
// The collision integral table is taken from Table E.2, and the viscosity
// equation is (1.4-14) of [1].
//
// Note that the conversion constant is calculated as
// 5/16*sqrt(k/N_a/pi) = 2.6695694574616934e-21 [in g/cm/s units]
// is slightly off the reported value in [1].
//
// Na = 6.02214129e23 molecules/mol
// k  =  1.3806488e-16 erg/K   (1.3806488e-23 J/K)
//
// [1] R.B. Bird, W.E. Stewart and E.N. Lightfood, Transport Phenomena,
//     2nd edition, J. Wiley, 2002. (accessed electronically)
const int KNumOmegaPoints = 82; 
const double KOmega2_2_LJ[] = {
  2.84,
  2.676,
  2.531,
  2.401,
  2.284,
  2.178,
  2.084,
  1.999,
  1.922,
  1.853,
  1.79,
  1.734,
  1.682,
  1.636,
  1.593,
  1.554,
  1.518,
  1.485,
  1.455,
  1.427,
  1.401,
  1.377,
  1.355,
  1.334,
  1.315,
  1.297,
  1.28,
  1.264,
  1.249,
  1.235,
  1.222,
  1.209,
  1.198,
  1.186,
  1.176,
  1.156,
  1.138,
  1.122,
  1.107,
  1.0933,
  1.0807,
  1.0691,
  1.0583,
  1.0482,
  1.0388,
  1.03,
  1.0217,
  1.0139,
  1.0066,
  0.9996,
  0.9931,
  0.9868,
  0.9809,
  0.9753,
  0.9699,
  0.9647,
  0.9598,
  0.9551,
  0.9506,
  0.9462,
  0.942,
  0.938,
  0.9341,
  0.9304,
  0.9268,
  0.8962,
  0.8727,
  0.8538,
  0.838,
  0.8244,
  0.8018,
  0.7836,
  0.7683,
  0.7552,
  0.7436,
  0.7198,
  0.701,
  0.6854,
  0.6723,
  0.651,
  0.614, 
  0.5887
};

const double KReducedTemperature[] = {
  0.3,
  0.35,
  0.4,
  0.45,
  0.5,
  0.55,
  0.6,
  0.65,
  0.7,
  0.75,
  0.8,
  0.85,
  0.9,
  0.95,
  1,
  1.05,
  1.1,
  1.15,
  1.2,
  1.25,
  1.3,
  1.35,
  1.4,
  1.45,
  1.5,
  1.55,
  1.6,
  1.65,
  1.7,
  1.75,
  1.8,
  1.85,
  1.9,
  1.95,
  2,
  2.1,
  2.2,
  2.3,
  2.4,
  2.5,
  2.6,
  2.7,
  2.8,
  2.9,
  3,
  3.1,
  3.2,
  3.3,
  3.4,
  3.5,
  3.6,
  3.7,
  3.8,
  3.9,
  4,
  4.1,
  4.2,
  4.3,
  4.4,
  4.5,
  4.6,
  4.7,
  4.8,
  4.9,
  5,
  6,
  7,
  8,
  9,
  10,
  12,
  14,
  16,
  18,
  20,
  25,
  30,
  35,
  40,
  50,
  75,
  100
};

} // end of namespace transport
#endif

