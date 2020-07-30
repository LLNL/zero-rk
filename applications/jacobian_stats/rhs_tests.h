#ifndef RHS_TESTS_H_
#define RHS_TESTS_H_

#include <zerork/mechanism.h>

int TestRateOfProgress(zerork::mechanism *mech,
                       const double temperature,
                       const double pressure,
                       const double mole_fraction[]);

#endif
