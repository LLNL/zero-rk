#ifndef ZERORK_CONSTANTS_H
#define ZERORK_CONSTANTS_H

#include <math.h>

namespace zerork {

const int MAX_ELEMENTS = 92;

//const double NIST_RU = 8.31447215e3;  // [J/kmol-k] gas constant in cantera
const double NIST_RU = 8.3144621e3;  // [J/kmol-k] accessed NIST site 7-8-2011
// http://physics.nist.gov/cgi-bin/cuu/Value?r
const double P_ATM = 1.01325e5;  // [Pa] reference pressure for Keq calculation



// multiplier to convert activation energies from cal/mol to kelvin
//const double CAL_PER_MOL_TACT = 4184.0/NIST_RU;
const double CAL_PER_MOL_TACT = 0.50321956485916268714484849236368519858909;


const int LDA_THERMO_POLY_D5R2 = 16; // minimum is 15 w/o (Tmin,Tmax) storage
const int NUM_THERMO_POLY_D5R2 = 15; // minimum is 15 w/o (Tmin,Tmax) storage
const int NUM_COEF_RANGE = 7; // number of coefficients per range

const int MIN_INT32 = -2147483648; // minimum 32-bit integer

const double STOICH_TOL = 1.0e-3; // threshold below which any stoichiometric
                                  // coefficient will be rounded to the nearest
                                  // integer and the reaction treated as an
                                  // elementary reaction

const double ARRHENIUS_RTOL = 1.0e-6;
const double ARRHENIUS_ATOL = 0.0;

const double MIN_EXP_ARG = log(1.0e-300);
const double MAX_EXP_ARG = log(1.0e+300);

const double KAvogadroNumber = 6.02214129e26; // [molecules/kmol]
const double KBoltzmann = 1.3806488e-23;      // [J/K]

 const int MAX_FILE_LINE = 1024;

//TODO: This is ugly.  Also need to have LA/HA
#ifdef ZERORK_USE_MKL
const char expType[]="intel_mkl";
#endif
#ifdef ZERORK_USE_FMATH
const char expType[]="fmath";
#endif
#ifdef ZERORK_USE_FMATH_NOCHECK
const char expType[]="fmath_nocheck";
#endif
#ifdef ZERORK_USE_LIBC
const char expType[]="libc";
#endif

} // namespace zerork

#endif
