#include <stdlib.h>
#include <stdio.h>

#include "special_species.h"

// getAdductSpecies returns the number of pure unimolecular species or adducts
// and a list of species indexes in the species_id vector, which will be
// cleared.  The current definition does not consider a species to be purely
// unimolecular if it is a reactant with a third body or falloff species or
// falloff third body. This definition may be modified in the future.
int getAdductSpeciesIndex(zerork::mechanism *mech,
                          std::vector<int> *species_id)
{
  const int num_species = mech->getNumSpecies();
  const int num_steps   = mech->getNumSteps();
  std::vector<int> max_reactant_order;
 
  max_reactant_order.assign(num_species,0);
  for(int j=0; j<num_steps; ++j) {
    int num_reactants = mech->getOrderOfStep(j);
    int reaction_id = mech->getRxnIdxOfStep(j);

    // determine the true reactant order of the step. Note that 
    // mech->getOrderOfStep() returns the number of reactant species ignoring
    // third body or falloff reactions.
    int reactant_order = num_reactants;
    if(mech->isThirdBodyReaction(reaction_id)) {
      reactant_order++;
    } else if(mech->isFalloffReaction(reaction_id)) {
      reactant_order++;
    }
    
    // set the maximum reactant order of a species
    for(int k=0; k<num_reactants; ++k) {
      int reactant_id = mech->getSpecIdxOfStepReactant(j,k);
      if(max_reactant_order[reactant_id] < reactant_order) {
        max_reactant_order[reactant_id] = reactant_order;
      }
    }
  }
  // At this point any pure unimolecular species will have a maximum
  // reactant order of one.  Non-reacting species will have a maximum
  // reactant order of zero.
  species_id->clear();
  for(int j=0; j<num_species; ++j) {
    if(max_reactant_order[j] == 1) {
      species_id->push_back(j);
    }
  }
  return static_cast<int>(species_id->size());
}

// getSourceSpecies returns the number of species that are "pure sources."
// That is, they are only reactants and never appear as products. In addition,
// a list of species indexes in the species_id vector will be updated as a
// a function argument. The search for pure sources does not check if a
// species is a special fall-off third-body species.  Note the species_id
// vector will be cleared, so that it will be empty if no species are found.
int getSourceSpeciesIndex(zerork::mechanism *mech,
                          std::vector<int> *species_id)
{
  const int num_species = mech->getNumSpecies();
  int num_source_species = 0;
  std::vector<int> reactant_count;
  std::vector<int> product_count;

  getSpeciesReactantCount(mech,&reactant_count);
  getSpeciesProductCount(mech,&product_count);

  for(int j=0; j<num_species; ++j) {
    if(reactant_count[j] > 0 && product_count[j] == 0) {
      species_id->push_back(j);
      ++num_source_species;
    }
  }

  return num_source_species;
}

// getSinkSpecies returns the number of species that are "pure sinks."
// That is, they are only products and never appear as reactants. In addition,
// a list of species indexes in the species_id vector will be updated as a
// a function argument. The search for pure sinks does not check if a
// species is a special fall-off third-body species.  Note the species_id
// vector will be cleared, so that it will be empty if no species are found.
int getSinkSpeciesIndex(zerork::mechanism *mech,
                        std::vector<int> *species_id)
{
  const int num_species = mech->getNumSpecies();
  int num_sink_species = 0;
  std::vector<int> reactant_count;
  std::vector<int> product_count;

  getSpeciesReactantCount(mech,&reactant_count);
  getSpeciesProductCount(mech,&product_count);

  for(int j=0; j<num_species; ++j) {
    if(reactant_count[j] == 0 && product_count[j] > 0) {
      species_id->push_back(j);
      ++num_sink_species;
    }
  }

  return num_sink_species;
}

int getSpeciesReactantCount(zerork::mechanism *mech,
                            std::vector<int> *reactant_count)
{
  const int num_species = mech->getNumSpecies();
  const int num_steps   = mech->getNumSteps();
  int species_id;
  int reactant_count_sum = 0;

  reactant_count->clear();
  reactant_count->assign(num_species,0);

  for(int j=0; j<num_steps; ++j) {

    int num_reactants = mech->getOrderOfStep(j);

    for(int k=0; k<num_reactants; ++k) {
      species_id = mech->getSpecIdxOfStepReactant(j,k);
      ++(*reactant_count)[species_id];
      ++reactant_count_sum;
    }
  }  
  return reactant_count_sum;
}

int getSpeciesProductCount(zerork::mechanism *mech,
                           std::vector<int> *product_count)
{
  const int num_species = mech->getNumSpecies();
  const int num_steps   = mech->getNumSteps();
  int species_id;
  int product_count_sum = 0;

  product_count->clear();
  product_count->assign(num_species,0);

  for(int j=0; j<num_steps; ++j) {

    int num_products = mech->getNumProductsOfStep(j);

    for(int k=0; k<num_products; ++k) {
      species_id = mech->getSpecIdxOfStepProduct(j,k);
      ++(*product_count)[species_id];
      ++product_count_sum;
    }
  }  
  return product_count_sum;
}

double getAdductRate(zerork::mechanism *mech, const int adduct_id,
                     const double temperature, const double pressure)
{
  const int num_species = mech->getNumSpecies();
  const int num_steps   = mech->getNumSteps();
  double adduct_concentration = pressure/(mech->getGasConstant()*temperature);
  std::vector<double> concentration;
  std::vector<double> net_reaction_rate;
  std::vector<double> creation_rate;
  std::vector<double> destruction_rate;
  std::vector<double> step_rate;

  concentration.assign(num_species,0.0);
  net_reaction_rate.assign(num_species,0.0);
  creation_rate.assign(num_species,0.0);
  destruction_rate.assign(num_species,0.0);
  step_rate.assign(num_steps,0.0);

  concentration[adduct_id] = adduct_concentration;

  mech->getReactionRates(temperature,
                         &concentration[0],
                         &net_reaction_rate[0],
                         &creation_rate[0],
                         &destruction_rate[0],
                         &step_rate[0]);

  return destruction_rate[adduct_id]/adduct_concentration;
}
