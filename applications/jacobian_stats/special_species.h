#ifndef SPECIAL_SPECIES_H_
#define SPECIAL_SPECIES_H_

#include "zerork/mechanism.h"

int getAdductSpeciesIndex(zerork::mechanism *mech,
                          std::vector<int> *species_id);
int getSourceSpeciesIndex(zerork::mechanism *mech,
                          std::vector<int> *species_id);
int getSinkSpeciesIndex(zerork::mechanism *mech,
                        std::vector<int> *species_id);

int getSpeciesReactantCount(zerork::mechanism *mech,
                            std::vector<int> *reactant_count);
int getSpeciesProductCount(zerork::mechanism *mech,
                           std::vector<int> *product_count);
double getAdductRate(zerork::mechanism *mech, const int adduct_id,
                     const double temperature, const double pressure);

#endif

