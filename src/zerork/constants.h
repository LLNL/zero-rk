#ifndef ZERORK_CONSTANTS_H
#define ZERORK_CONSTANTS_H

#include <math.h>

namespace zerork {

const int MAX_ELEMENTS = 92;

const double NIST_RU = 8.314462618e3; // [J/kmol-k]; accessed NIST site 11-22-2024

const double P_ATM = 1.01325e5;  // [Pa] reference pressure for Keq calculation


const double KAvogadroNumber = 6.02214076e26; // [molecules/kmol]
const double KBoltzmann = 1.380649e-23; // [J/K]

// multiplier to convert activation energies from cal/mol to kelvin
const double CAL_PER_MOL_TACT = 4184.0/NIST_RU;
const double KCAL_PER_MOL_TACT = 4184.0/NIST_RU*1000.0;

// convert from kj/mol to cal/mol then cal/mol to kelvin
const double KJOULES_PER_MOL_TACT = 238.846*CAL_PER_MOL_TACT;
const double JOULES_PER_MOL_TACT = 238.846*CAL_PER_MOL_TACT/1000.0;

const int LDA_THERMO_POLY_D5R2 = 16; // minimum is 15 w/o (Tmin,Tmax) storage
const int NUM_THERMO_POLY_D5R2 = 15; // minimum is 15 w/o (Tmin,Tmax) storage
const int NUM_COEF_RANGE = 7; // number of coefficients per range

const int MIN_INT32 = (-2147483647-1); // minimum 32-bit integer

const double STOICH_TOL = 1.0e-14; // threshold below which any stoichiometric
                                   // coefficient will be rounded to the nearest
                                   // integer and the reaction treated as an
                                   // elementary reaction

const double ARRHENIUS_RTOL = 1.0e-6;
const double ARRHENIUS_ATOL = 0.0;

const double MIN_EXP_ARG = log(1.0e-300);
const double MAX_EXP_ARG = log(1.0e+300);

const int MAX_FILE_LINE = 1024;

} // namespace zerork

#endif
