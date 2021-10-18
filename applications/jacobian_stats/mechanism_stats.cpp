#include <stdlib.h>
#include <stdio.h>
#include <limits>

#include "utilities/sort_vector.h"

#include "special_species.h"
#include "mechanism_stats.h"


int getBasicMechanismReport(zerork::mechanism *mech, std::string *report)
{
  char format_line[MAX_LINE_LENGTH];
  std::vector<int> step_type_count;

  report->clear();
  *report = "# Basic Mechanism Report\n";

  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of species       : %d\n",
           mech->getNumSpecies());  
  (*report)+=std::string(format_line);

  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of reactions     : %d\n",
           mech->getNumReactions());  
  (*report)+=std::string(format_line);

  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of 1-way steps   : %d\n",
           mech->getNumSteps());  
  (*report)+=std::string(format_line);

  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   max number of reactants : %d\n",
           mech->getMaxReactantsInStep());  
  (*report)+=std::string(format_line);

  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   max number of products  : %d\n",
           mech->getMaxProductsInStep());  
  (*report)+=std::string(format_line);

  getStepTypeCount(mech,&step_type_count);
  for(unsigned int j=0; j<step_type_count.size(); ++j) {

    int num_reactants = (j/mech->getMaxProductsInStep())+1;
    int num_products  = (j%mech->getMaxProductsInStep())+1;

     snprintf(format_line,
              MAX_LINE_LENGTH,
              "#     number of (%d) => (%d) steps : %d\n",
              num_reactants,num_products,step_type_count[j]);
     (*report)+=std::string(format_line);
  }

  return 0;
}

int getSpeciesSummaryReport(zerork::mechanism *mech, std::string *report)
{
  char format_line[MAX_LINE_LENGTH];
  int species_count;
  std::vector<int> species_id;

  report->clear();
  *report = "# Species Summary Report:\n";

  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   total number of species         : %d\n",
           mech->getNumSpecies());  
  (*report)+=std::string(format_line);

  species_count = getInertSpeciesIndex(mech,&species_id);

  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of inert species         : %d\n",
           species_count);  
  (*report)+=std::string(format_line);

  species_count = getLoneReactionSpeciesIndex(mech,&species_id);

  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of lone reaction species : %d\n",
           species_count);  
  (*report)+=std::string(format_line);

  species_count = getAdductSpeciesIndex(mech,&species_id);

  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of adduct species        : %d\n",
           species_count);  
  (*report)+=std::string(format_line);

  species_count = getSourceSpeciesIndex(mech,&species_id);

  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of source species        : %d\n",
           species_count);  
  (*report)+=std::string(format_line);

  species_count = getSinkSpeciesIndex(mech,&species_id);

  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of sink species          : %d\n",
           species_count);  
  (*report)+=std::string(format_line);


  return static_cast<int>(report->size());
}

int getInertSpeciesIndex(zerork::mechanism *mech,
                          std::vector<int> *species_id)
{
  int num_inert_species=0;
  std::vector<int> step_count;

  //step_count.resize(mech->getNumSpecies());

  buildSpeciesStepCount(mech,
                        &step_count);
  species_id->clear();

  for(int j=0; j<mech->getNumSpecies(); ++j) {
    if(step_count[j] == 0) { // found an inert species
      species_id->push_back(j);
      ++num_inert_species;
    }
  }
  return num_inert_species;
}

void getInertSpeciesReport(zerork::mechanism *mech,
                           std::string *report)
{
  int num_inert;
  std::vector<int> species_id;
  char format_line[MAX_LINE_LENGTH];

  num_inert=getInertSpeciesIndex(mech,
                                 &species_id);
  report->clear();
  *report = "# Inert Species Report\n";
 
  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of inert species : %d\n",
           num_inert);  
  (*report)+=std::string(format_line);
  if(num_inert > 0) {
    (*report)+="#   index, name (species index starts at 0)\n";
  }
  for(int j=0; j<num_inert; ++j) {
    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#  %6d  %s\n",
             species_id[j],
             mech->getSpeciesName(species_id[j]));
    (*report)+=std::string(format_line);
  }   
}

int getLoneReactionSpeciesIndex(zerork::mechanism *mech,
                                std::vector<int> *species_id)
{
  int num_lone_species=0;
  std::vector<int> reaction_count;

  buildSpeciesReactionCount(mech,
                            &reaction_count);
  species_id->clear();
  for(int j=0; j<mech->getNumSpecies(); ++j) {
    if(reaction_count[j] == 1) { // found a species in only one reaction
      species_id->push_back(j);
      
      ++num_lone_species;
    }
  }
  return num_lone_species;
}

void getLoneReactionSpeciesReport(zerork::mechanism *mech,
                                  std::string *report)
{
  int num_lone;
  std::vector<int> species_id;
  std::vector<int> reaction_id;
  char format_line[MAX_LINE_LENGTH];

  num_lone=getLoneReactionSpeciesIndex(mech,
                                       &species_id);
  report->clear();
  *report = "# Lone Reaction Species Report\n";
 
  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of species in only one reaction : %d\n",
           num_lone);  
  (*report)+=std::string(format_line);
  if(num_lone > 0) {
    (*report)+="#   Note that the species/reaction index starts at 0.\n";
    (*report)+="#   index, name, reaction index, reaction\n";
  }
  for(int j=0; j<num_lone; ++j) {

    buildReactionListOfSpecies(species_id[j],
                               mech,
                               &reaction_id);

    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#  %6d  %16s  %6d  %s\n",
             species_id[j],
             mech->getSpeciesName(species_id[j]),
             reaction_id[0],
             mech->getReactionString(reaction_id[0]).c_str());

    (*report)+=std::string(format_line);
  }  
}

void getAdductSpeciesReport(zerork::mechanism *mech,
                           std::string *report)
{
  int num_adduct;
  std::vector<int> species_id;
  char format_line[MAX_LINE_LENGTH];

  num_adduct=getAdductSpeciesIndex(mech,
                                   &species_id);
  report->clear();
  *report  = "# Adduct Species Report\n";
  *report += "# An \"adduct species\" is defined as a species that is the only reactant\n";
  *report += "# in all reaction steps in which it is involved. Here this excludes\n";
  *report += "# reaction steps with a single reactant interacting with a third-body\n";
  *report += "# falloff or otherwise.\n"; 
  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of adduct species : %d\n",
           num_adduct);  
  (*report)+=std::string(format_line);
  if(num_adduct > 0) {
    (*report)+="#   species index starts at 0| adduct destruction rate [Hz] at given temps  |\n";
    (*report)+="#   index, name              300 K     600 K     1200 K    1800 K    3000 K\n";
  }
  for(int j=0; j<num_adduct; ++j) {
    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#  %6d  %-16s  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e\n",
             species_id[j],
             mech->getSpeciesName(species_id[j]),
             getAdductRate(mech,species_id[j], 300.0,1.01325e5),
             getAdductRate(mech,species_id[j], 600.0,1.01325e5),
             getAdductRate(mech,species_id[j],1200.0,1.01325e5),
             getAdductRate(mech,species_id[j],1800.0,1.01325e5),
             getAdductRate(mech,species_id[j],3000.0,1.01325e5));

    (*report)+=std::string(format_line);
  }   
}

void getSourceSpeciesReport(zerork::mechanism *mech,
                            std::string *report)
{
  int num_source;
  std::vector<int> species_id;
  char format_line[MAX_LINE_LENGTH];

  num_source=getSourceSpeciesIndex(mech,
                                   &species_id);
  report->clear();
  *report  = "# Source Species Report\n";
  *report += "# A \"source species\" is defined as a species that only appears as a\n";
  *report += "# reactant. It may be defined in the initial composition, but will only\n";
  *report += "# decrease in quantity. It will never be produced.\n";
  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of source species : %d\n",
           num_source);  
  (*report)+=std::string(format_line);
  if(num_source > 0) {
    (*report)+="#   species index starts at 0\n";
    (*report)+="#   index, name\n";
  }
  for(int j=0; j<num_source; ++j) {
    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#  %6d  %-16s\n",
             species_id[j],
             mech->getSpeciesName(species_id[j]));

    (*report)+=std::string(format_line);
  }   
}

void getSinkSpeciesReport(zerork::mechanism *mech,
                          std::string *report)
{
  int num_sink;
  std::vector<int> species_id;
  char format_line[MAX_LINE_LENGTH];

  num_sink=getSinkSpeciesIndex(mech,
                                   &species_id);
  report->clear();
  *report  = "# Sink Species Report\n";
  *report += "# A \"sink species\" is defined as a species that only appears as a\n";
  *report += "# product. It  will only increase in quantity, and will never be consumed.\n";
  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of sink species : %d\n",
           num_sink);  
  (*report)+=std::string(format_line);
  if(num_sink > 0) {
    (*report)+="#   species index starts at 0\n";
    (*report)+="#   index, name\n";
  }
  for(int j=0; j<num_sink; ++j) {
    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#  %6d  %-16s\n",
             species_id[j],
             mech->getSpeciesName(species_id[j]));

    (*report)+=std::string(format_line);
  }   
}

void buildSpeciesReactionCount(zerork::mechanism *mech,
                        std::vector<int> *reaction_count)
{
  reaction_count->assign(mech->getNumSpecies(),0);

  for(int j=0; j<mech->getNumReactions(); ++j) {

    int step_id = mech->getStepIdxOfRxn(j,1);

    for(int k=0; k<mech->getOrderOfStep(step_id); ++k) {
      ++(*reaction_count)[mech->getSpecIdxOfStepReactant(step_id,k)];
    }

    for(int k=0; k<mech->getNumProductsOfStep(step_id); ++k) {
      ++(*reaction_count)[mech->getSpecIdxOfStepProduct(step_id,k)];
    }

  }
}

void buildSpeciesStepCount(zerork::mechanism *mech,
                           std::vector<int> *step_count)
{
  step_count->assign(mech->getNumSpecies(),0);

  for(int j=0; j<mech->getNumSteps(); ++j) {
    for(int k=0; k<mech->getOrderOfStep(j); ++k) {
      ++(*step_count)[mech->getSpecIdxOfStepReactant(j,k)];
    }
    for(int k=0; k<mech->getNumProductsOfStep(j); ++k) {
      ++(*step_count)[mech->getSpecIdxOfStepProduct(j,k)];
    }
  }
  
}

int buildReactionListOfSpecies(const int species_id,
                               zerork::mechanism *mech, 
                               std::vector<int> *reaction_list)
{
  int num_found=0;
  reaction_list->clear();

  for(int j=0; j<mech->getNumReactions(); ++j) {

    int step_id = mech->getStepIdxOfRxn(j,1);
    int is_target = 0;

    for(int k=0; k<mech->getOrderOfStep(step_id); ++k) {
      if(mech->getSpecIdxOfStepReactant(step_id,k) == species_id) {
        is_target=1;
      }
    }

    for(int k=0; k<mech->getNumProductsOfStep(step_id); ++k) {
      if(mech->getSpecIdxOfStepProduct(step_id,k) == species_id) {
        is_target=1;
      }
    }
    if(is_target==1) {
      reaction_list->push_back(j);
      ++num_found;
    }
  }
    
  return num_found;
}


int updateSpeciesActivityOfStep(const int step_id,
                                zerork::mechanism *mech,
                                std::vector<int> *active_species)
{
  int is_step_active=0;

  if(static_cast<int>(active_species->size()) != mech->getNumSpecies()) {
    printf("ERROR: In updateSpeciesActivityOfStep(...),\n");
    printf("       active_species length (%lu) does not equal\n",
           active_species->size());
    printf("       the total species in the mechanism (%d)\n",
           mech->getNumSpecies());
    exit(-1);
  }
  // if first species is present (active) set the active flag
  if((*active_species)[mech->getSpecIdxOfStepReactant(step_id,0)] > 0) {
    is_step_active = 1;
  }
  // check the remaining species presence
  for(int k=1; k<mech->getOrderOfStep(step_id); ++k) {
    is_step_active *=
      (*active_species)[mech->getSpecIdxOfStepReactant(step_id,k)];
  }
  if(is_step_active > 0) { // activate the products
    for(int k=0; k<mech->getNumProductsOfStep(step_id); ++k) {
      (*active_species)[mech->getSpecIdxOfStepProduct(step_id,k)]=1;
    }
  }
  return is_step_active;
}

int getInactiveSpeciesAndSteps(std::vector<int> *initial_active_species,
                               zerork::mechanism *mech,
                               std::vector<int> *inactive_species,
                               std::vector<int> *inactive_steps)
{
  int num_active_species;
  int prev_active_steps=-1;
  int num_active_steps;
  int num_generations=0;
  std::vector<int> active_species;
  std::vector<int> active_steps;
  
  active_species.assign(mech->getNumSpecies(),0);
  for(unsigned int j=0; j<initial_active_species->size(); ++j) {
    active_species[(*initial_active_species)[j]] = 1;
  }

  num_active_species=initial_active_species->size();
  num_active_steps=0;
  active_steps.assign(mech->getNumSteps(),0);

  while(num_active_steps != prev_active_steps) {

    prev_active_steps = num_active_steps;
    num_active_steps = 0;
    for(int j=0; j<mech->getNumSteps(); ++j) {
 
      active_steps[j] = updateSpeciesActivityOfStep(j,
                                                    mech,
                                                    &active_species);
      num_active_steps += active_steps[j];
    }

    // count the number of active species at this generation
    num_active_species=0;
    for(int j=0; j<mech->getNumSpecies(); ++j) {
      num_active_species+=active_species[j];
    }

    ++num_generations;  
    //printf("generation %d: active species %d, active steps %d\n",
    //       num_generations,num_active_species, num_active_steps);
  }
  --num_generations; // last generation didn't produce any newly active steps

  inactive_species->clear();
  for(int j=0; j<mech->getNumSpecies(); ++j) {
    if(active_species[j] == 0) {
      inactive_species->push_back(j);
    }
  }

  inactive_steps->clear();
  for(int j=0; j<mech->getNumSteps(); ++j) {
    if(active_steps[j] == 0) {
      inactive_steps->push_back(j);
    }
  }

  return num_generations;
}

void getInactiveReport(std::vector<int> *initial_active_species,
                       zerork::mechanism *mech,
                       std::string *report)
{
  int num_generations;
  std::vector<int> inactive_species_id;
  std::vector<int> inactive_steps_id;
  char format_line[MAX_LINE_LENGTH];

  num_generations=getInactiveSpeciesAndSteps(initial_active_species,
                                             mech,
                                             &inactive_species_id,
                                             &inactive_steps_id);
  report->clear();

  (*report)  = "# Inactive Species and Steps Report\n";
  (*report) += "#   initial active species composition:\n";
  (*report) += "#     index, name\n";

  for(unsigned int j=0; j<initial_active_species->size(); ++j) {
    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#     %6d  %s\n",
             (*initial_active_species)[j],
	     mech->getSpeciesName((*initial_active_species)[j]));  
    (*report)+=std::string(format_line);
  }
 
  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   reaction network maximally active at generation %d\n",
           num_generations);
  (*report)+=std::string(format_line);  
 
  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of inactive species = %lu\n",
           inactive_species_id.size());
  (*report)+=std::string(format_line);  

  if(inactive_species_id.size()>0) {
    (*report) += "#     index, name\n";
  }

  for(unsigned int j=0; j<inactive_species_id.size(); ++j) {
    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#     %6d  %s\n",
             inactive_species_id[j],
             mech->getSpeciesName(inactive_species_id[j]));
    (*report)+=std::string(format_line);  
  }    
  
  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of inactive 1-way steps = %lu\n",
           inactive_steps_id.size());
  (*report)+=std::string(format_line);  

  if(inactive_steps_id.size()>0) {
    (*report) += "#     step index, reaction index, direction, name\n";
  }

  for(unsigned int j=0; j<inactive_steps_id.size(); ++j) {
    int reaction_id = mech->getRxnIdxOfStep(inactive_steps_id[j]);
    char step_direction[] = "(fwd)";
    if(mech->getStepIdxOfRxn(reaction_id,-1) == inactive_steps_id[j]) {
      strcpy(step_direction,"(rev)");
    }
    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#     %6d  %6d  %s  %s\n",
             inactive_steps_id[j],
             reaction_id,
             step_direction,
             mech->getReactionString(reaction_id).c_str());
    (*report)+=std::string(format_line);  
  }
}

// The step_type_count vector is ordered so that the product counts increment
// in the inner-loop.  That is,
//
// step_type_count[0] = count of 1 reactant -> 1 product steps
// step_type_count[1] = count of 1 reactant -> 2 product steps
//                     :
// step_type_count[max_prod] = count of 2 reactant -> 1 product steps
//                            :
// where max_prod is the maximum number of products in all steps
void getStepTypeCount(zerork::mechanism *mech,
                       std::vector<int> *step_type_count)
{
  step_type_count->clear();
  step_type_count->assign(mech->getMaxReactantsInStep()*
                          mech->getMaxProductsInStep(),
			  0);
  for(int j=0; j<mech->getNumSteps(); ++j) {
    int order = mech->getOrderOfStep(j);
    int prods = mech->getNumProductsOfStep(j);
    if(order <= 0 || prods <= 0) {
      //Non-integer step
      continue;
    }
    int bin_id = (order-1)*mech->getMaxProductsInStep()+(prods-1);
    ++(*step_type_count)[bin_id];
  }
}

void scanNegativeRate(zerork::mechanism *mech,
                      const double pressure_ref,
                      const double min_temperature,
                      const double max_temperature,
                      const int num_scans,
                      int *min_rate_step_id,
                      double *min_rate_temperature,
                      double *min_rate)
{
  double current_temperature;
  double average_concentration;
  std::vector<int> forward_id;
  std::vector<int> reverse_id;
  *min_rate = std::numeric_limits<double>::max();

  // build the list of forward and reverse reaction indexes
  for(int j=0; j<mech->getNumSteps(); ++j) {
    int reaction_id = mech->getRxnIdxOfStep(j);
    if(mech->getStepIdxOfRxn(reaction_id,1) == j) {
      forward_id.push_back(reaction_id);
    }
    else {
      reverse_id.push_back(reaction_id);
    }
  }
  
  // Create an array of dummy concentrations to calculate the forward and
  // reverse rate constants.  Since the third body and falloff reactions are
  // not considered unimolecular, the dummy concentrations should not affect
  // the search for the fastest unimolecular rates.
  std::vector<double> dummy_concentration(mech->getNumSpecies());
  // Create arrays to store the forward and reverse rate constants
  std::vector<double> forward_rate_const(mech->getNumReactions());
  std::vector<double> reverse_rate_const(mech->getNumReactions());

  for(int j=0; j<num_scans; ++j) {

    current_temperature = min_temperature + (max_temperature-min_temperature)*
      static_cast<double>(j)/static_cast<double>(num_scans-1);

    average_concentration = pressure_ref/(mech->getGasConstant()*current_temperature)/
      static_cast<double>(mech->getNumSpecies());

    for(int j=0; j<mech->getNumSpecies(); ++j) {
      dummy_concentration[j] = average_concentration;
    }

    mech->getKrxnFromTC(current_temperature,
                        &dummy_concentration[0],
                        &forward_rate_const[0],
                        &reverse_rate_const[0]);
    // check the unimolecular rate constants in the forward direction
    for(unsigned int k=0; k<forward_id.size(); ++k) {
      if(forward_rate_const[forward_id[k]] < *min_rate) {
        (*min_rate_temperature) = current_temperature;
        (*min_rate_step_id) = mech->getStepIdxOfRxn(forward_id[k],-1);
        (*min_rate) = forward_rate_const[forward_id[k]];
      }
    }
    // check the unimolecular rate constants in the reverse direction
    for(unsigned int k=0; k<reverse_id.size(); ++k) {
      if(reverse_rate_const[reverse_id[k]] < *min_rate) {
        (*min_rate_temperature) = current_temperature;
        (*min_rate_step_id) = mech->getStepIdxOfRxn(reverse_id[k],-1);
        (*min_rate) = reverse_rate_const[reverse_id[k]];
      }
    }
  } // j-loop over the number of scans
}

double scanMaxUnimolecularRate(zerork::mechanism *mech,
                               const double pressure_ref,
                               const double min_temperature,
                               const double max_temperature,
                               const int num_scans,
                               double *max_rate_temperature,
                               int *max_rate_step_id)
{
  double current_temperature;
  double average_concentration;
  double *dummy_concentration;
  double *forward_rate_const,*reverse_rate_const;
  double max_unimolecular_rate=0.0;
  int num_unimolecular_steps;
  std::vector<int> forward_id;
  std::vector<int> reverse_id;
  std::vector<int> unimolecular_step_id;
  forward_id.clear();
  reverse_id.clear();

  num_unimolecular_steps = buildUnimolecularListOfSteps(mech,
                                                        &unimolecular_step_id);
  if(num_unimolecular_steps < 1) {
    printf("WARNING: In scanMaxUnimolecularRate(...),\n");
    printf("         the number of unimolecular steps is %d.\n",
           num_unimolecular_steps);
    max_rate_temperature = NULL;
    max_rate_step_id = NULL;
    return 0.0;
  }

  // build the list of forward and reverse reaction indexe corresponding to 
  // true unimolecular reactions steps                                         
  for(unsigned int j=0; j<unimolecular_step_id.size(); ++j) {
    int reaction_id = mech->getRxnIdxOfStep(unimolecular_step_id[j]);
    if(mech->getStepIdxOfRxn(reaction_id,1) == unimolecular_step_id[j]) {
      forward_id.push_back(reaction_id);
    }
    else {
      reverse_id.push_back(reaction_id);
    }
  }
  
  // Create an array of dummy concentrations to calculate the forward and
  // reverse rate constants.  Since the third body and falloff reactions are
  // not considered unimolecular, the dummy concentrations should not affect
  // the search for the fastest unimolecular rates.
  dummy_concentration = new double[mech->getNumSpecies()];
  // Create arrays to store the forward and reverse rate constants
  forward_rate_const = new double[mech->getNumReactions()];
  reverse_rate_const = new double[mech->getNumReactions()];


  for(int j=0; j<num_scans; ++j) {

    current_temperature = min_temperature + (max_temperature-min_temperature)*
      static_cast<double>(j)/static_cast<double>(num_scans-1);

    average_concentration = pressure_ref/(mech->getGasConstant()*current_temperature)/
      static_cast<double>(mech->getNumSpecies());
  
    for(int j=0; j<mech->getNumSpecies(); ++j) {
      dummy_concentration[j] = average_concentration;
    }
    mech->getKrxnFromTC(current_temperature,
                        &dummy_concentration[0],
                        &forward_rate_const[0],
                        &reverse_rate_const[0]);
    // check the unimolecular rate constants in the forward direction
    for(unsigned int k=0; k<forward_id.size(); ++k) {
      if(forward_rate_const[forward_id[k]] > max_unimolecular_rate) {
        (*max_rate_temperature) = current_temperature;
        (*max_rate_step_id) = mech->getStepIdxOfRxn(forward_id[k],1);
        max_unimolecular_rate = forward_rate_const[forward_id[k]];
      }
    }
    // check the unimolecular rate constants in the reverse direction
    for(unsigned int k=0; k<reverse_id.size(); ++k) {
      if(reverse_rate_const[reverse_id[k]] > max_unimolecular_rate) {
        (*max_rate_temperature) = current_temperature;
        (*max_rate_step_id) = mech->getStepIdxOfRxn(reverse_id[k],-1);
        max_unimolecular_rate = reverse_rate_const[reverse_id[k]];
      }
    }
  } // j-loop over the number of scans
    

  delete [] forward_rate_const;
  delete [] reverse_rate_const;
  delete [] dummy_concentration;

  return max_unimolecular_rate;
}

int buildUnimolecularListOfSteps(zerork::mechanism *mech,
				 std::vector<int> *step_id)
{
  step_id->clear();
  for(int j=0; j<mech->getNumSteps(); ++j) {
    int reaction_id = mech->getRxnIdxOfStep(j);
    if(mech->isThirdBodyReaction(reaction_id) == 0
       && mech->isFalloffReaction(reaction_id) == 0) {
      if(mech->getOrderOfStep(j) == 1) {
        step_id->push_back(j);
      }
    }
  }
  return step_id->size();
}

void getNegativeRateReport(zerork::mechanism *mech,
                           const double pressure_ref,
                           const double min_temperature,
                           const double max_temperature,
			   const int num_scans,
                           std::string *report)
{
  char format_line[MAX_LINE_LENGTH];
  char step_direction[] = "(fwd)";
  int min_rate_step_id;
  double min_rate_temperature;
  double min_rate;

  report->clear();
  (*report) = "# Negative Rate Report\n";

  if(mech->getNumSteps() > 0) {
    scanNegativeRate(mech,
                     pressure_ref,
                     min_temperature,
                     max_temperature,
                     num_scans,
                     &min_rate_step_id,
                     &min_rate_temperature,
                     &min_rate);

    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#   scanned temperature [%7.2f,%7.2f] in increments of %10f\n",
             min_temperature,
             max_temperature,
             (max_temperature-min_temperature)/static_cast<double>(num_scans));
    (*report)+=std::string(format_line);

    if(min_rate >= 0) {
      snprintf(format_line,
               MAX_LINE_LENGTH,
               "#     no negative rates found\n");  
      (*report)+=std::string(format_line);
    } else {
      snprintf(format_line,
               MAX_LINE_LENGTH,
               "#     most negative rate:  %14.7e\n",
               min_rate); 
      (*report)+=std::string(format_line);

      snprintf(format_line,
               MAX_LINE_LENGTH,
               "#     at temperature          [K]:  %10.4f\n",
               min_rate_temperature);  
      (*report)+=std::string(format_line);

      int reaction_id = mech->getRxnIdxOfStep(min_rate_step_id);

      if(mech->getStepIdxOfRxn(reaction_id,-1) == min_rate_step_id) {
        strcpy(step_direction,"(rev)");
      }

      snprintf(format_line,
               MAX_LINE_LENGTH,
               "#     for step %6d            : %s %s\n",
               min_rate_step_id,
               mech->getReactionString(reaction_id).c_str(),
               step_direction);  
      (*report)+=std::string(format_line);
    }
  }
}

void getMaxUnimolecularRateReport(zerork::mechanism *mech,
                                  const double pressure_ref,
                                  const double min_temperature,
                                  const double max_temperature,
				  const int num_scans,
                                  std::string *report)
{
  char format_line[MAX_LINE_LENGTH];
  std::vector<int> unimolecular_step_id;
  double max_rate;
  double max_rate_temperature;
  int max_rate_step_id;
  int reaction_id;
  char step_direction[] = "(fwd)";

  report->clear();
  (*report) = "# Maximum Unimolecular Rate Report\n";

  buildUnimolecularListOfSteps(mech,&unimolecular_step_id);
  snprintf(format_line,
           MAX_LINE_LENGTH,
           "#   number of unimolecular steps : %lu (excluding 3rd-body & falloff)\n",
           unimolecular_step_id.size());  
  (*report)+=std::string(format_line);

  if(unimolecular_step_id.size() > 0) {
    max_rate = scanMaxUnimolecularRate(mech,
                                       pressure_ref,
                                       min_temperature,
                                       max_temperature,
                                       num_scans,
                                       &max_rate_temperature,
                                       &max_rate_step_id);

    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#   scanned temperature [%7.2f,%7.2f] in increments of %10f\n",
             min_temperature,
             max_temperature,
             (max_temperature-min_temperature)/static_cast<double>(num_scans));     (*report)+=std::string(format_line);
    
    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#     max unimolecular rate [1/s]:  %14.7e\n",
             max_rate);  
    (*report)+=std::string(format_line);

    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#     at temperature          [K]:  %10.4f\n",
             max_rate_temperature);  
    (*report)+=std::string(format_line);

    reaction_id = mech->getRxnIdxOfStep(max_rate_step_id);
    
    if(mech->getStepIdxOfRxn(reaction_id,-1) == max_rate_step_id) {
      strcpy(step_direction,"(rev)");
    }

    snprintf(format_line,
             MAX_LINE_LENGTH,
             "#     for step %6d            : %s %s\n",
             max_rate_step_id,
	     mech->getReactionString(reaction_id).c_str(),
             step_direction);  
    (*report)+=std::string(format_line);
    
  }

}
void getSortedUnimolecularRateReport(zerork::mechanism *mech,
                                     const double pressure_ref,
                                     const double min_temperature,
                                     const double max_temperature,
                                     const int num_scans,
                                     const double min_report_rate,
                                     std::string *report)
{
  double dummy_concentration;
  int reaction_id;
  int num_uni_steps;
  char format_line[zerork::MAX_FILE_LINE];
  int num_reported_steps;
  std::string step_name;
  zerork::utilities::SortableVector initial;
  std::vector<int> uni_step_id; 
  std::vector<double> temperature;
  std::vector<double> concentration, Kfwd, Krev;
  std::vector<zerork::utilities::SortableVector> uni_rates;
 
  // initialization
  temperature.assign(num_scans,0.0);
  concentration.assign(mech->getNumSpecies(),dummy_concentration);
  Kfwd.assign(mech->getNumReactions(),0.0);
  Krev.assign(mech->getNumReactions(),0.0);

  for(int j=0; j<num_scans; j++) {
    temperature[j] = min_temperature + (max_temperature-min_temperature)*
      (static_cast<double>(j)/static_cast<double>(num_scans-1));
  }
  num_uni_steps=buildUnimolecularListOfSteps(mech,&uni_step_id);
  if(num_uni_steps==0) return;

  for(int j=0; j<num_uni_steps; ++j) {
    initial.id = uni_step_id[j];
    initial.v_double.assign(num_scans,0.0);
    uni_rates.push_back(initial); 
  }
  // compute each unimolecular rate
  for(int j=0; j<num_scans; j++) {
  
    dummy_concentration = pressure_ref/
      (mech->getGasConstant()*min_temperature)/
      static_cast<double>(mech->getNumSpecies());
    concentration.assign(mech->getNumSpecies(),
                         dummy_concentration);
    mech->getKrxnFromTC(temperature[j],
                        concentration.data(),
                        Kfwd.data(),
                        Krev.data());
    for(int k=0; k<num_uni_steps; k++) {
      reaction_id = mech->getRxnIdxOfStep(uni_rates[k].id);
      if(mech->getStepIdxOfRxn(reaction_id,1) == uni_rates[k].id) {
        // forward rate
        uni_rates[k].v_double[j] = Kfwd[reaction_id];
      } else {
        uni_rates[k].v_double[j] = Krev[reaction_id];
      }
    }        
  }
  SortByVectorMax(num_uni_steps,uni_rates.data());
  
  // using the sorted rate data to determine how many steps exceed the
  // user threshold
  num_reported_steps = 0;
  for(int j=0; j<num_uni_steps; ++j) {
    if(SortableVectorMax(&uni_rates[j]) < min_report_rate) {
      num_reported_steps = j;
      break;
    }
  }

  report->clear();
  *report = "# Sorted Unimolecular Reaction Rate Report:\n";
  if(num_scans >= 2) {
    // range scan
    snprintf(format_line,zerork::MAX_FILE_LINE,
             "#   Highest rate found over T-range [%7.2f K,%7.2f K]\n",
             min_temperature,
             max_temperature);
    *report += std::string(format_line);
    snprintf(format_line,zerork::MAX_FILE_LINE,
             "#   scanning temperature steps dT = %9.4f K\n",
             temperature[1]-temperature[0]);
    *report += std::string(format_line);
  } else {
    // single point 
     snprintf(format_line,zerork::MAX_FILE_LINE,
             "#   Highest rate found for T = %7.2f K\n",
	      temperature[0]);   
    *report += std::string(format_line);
  }
  snprintf(format_line,zerork::MAX_FILE_LINE,
           "#   and p = %12.4e Pa (will only impact PLOG reactions)\n",
	    pressure_ref);
  *report += std::string(format_line);

  mech->getReactionNameDirOfStep(uni_rates[0].id, &step_name);

  snprintf(format_line,zerork::MAX_FILE_LINE,
           "#   fastest unimolecular step: %d {%s}\n",
           uni_rates[0].id,
           step_name.c_str());
  *report += std::string(format_line);

  snprintf(format_line,zerork::MAX_FILE_LINE,
           "#   step %7d rate [1/s]: %12.5e\n",
           uni_rates[0].id,
           SortableVectorMax(&uni_rates[0]));
  *report += std::string(format_line);

  snprintf(format_line,zerork::MAX_FILE_LINE,
           "#   Number of unimolecur steps reported with a peak rate\n");
  *report += std::string(format_line);

  snprintf(format_line,zerork::MAX_FILE_LINE,
           "#   greater than or equal to %8.2e [1/s]: %d\n",
           min_report_rate,
           num_reported_steps);
  *report += std::string(format_line);

  if(num_reported_steps < 1) {
    return;
  }

  if(num_scans == 1) {
    // write out the single temperature probability as a list
    for(int j=0; j<num_reported_steps; ++j) {
      mech->getReactionNameDirOfStep(uni_rates[j].id, &step_name);
      snprintf(format_line,zerork::MAX_FILE_LINE,
               "%5d  %10.3e { %s }\n",
               uni_rates[j].id,
               uni_rates[j].v_double[0],
               step_name.c_str());
      *report += std::string(format_line); 
    }
  } else {
    // write out the
    // TODO: apply threshold to reduce the number of unimolecular rate reports
    WriteUnimolecularScanInBlocks(mech,
                                  7,
                                  num_scans,
                                  num_reported_steps,
                                  temperature.data(),
                                  &uni_rates,
                                  report);
  } 
}
void WriteUnimolecularScanInBlocks(zerork::mechanism *mech,
                           const int block_size,
                           const int num_scans,
			   const int num_uni_steps,
                           const double scan_temperature[],
			   const std::vector<zerork::utilities::SortableVector> *probability,
			   std::string *report)
{
  int step_id, list_id;
  char format_line[zerork::MAX_FILE_LINE];
  std::string step_name;
  int num_blocks  = num_uni_steps/block_size;
  int rem_columns = num_uni_steps%block_size;

  for(int j=0; j<num_blocks; ++j) {
    // write new block header
    snprintf(format_line,zerork::MAX_FILE_LINE,
             "# gnuplot data block (index): %d\n",
             j);
    *report += std::string(format_line);
    *report += "# column 1: [K] scan temperature\n";
    for(int k=0; k<block_size; ++k) {
      list_id = k+block_size*j;
      step_id = probability->at(list_id).id;
      mech->getReactionNameDirOfStep(step_id, &step_name);

      snprintf(format_line,zerork::MAX_FILE_LINE,
               "# column %d: [1/s] step %4d rate {max =%10.3e %s }\n",
               k+2,
               step_id,
               SortableVectorMax(&probability->at(list_id)),
               step_name.c_str());
      *report += std::string(format_line);
    }
    *report += "#------------------------------------------------------------------------------\n";
    // end new block header
    for(int k=0; k<num_scans; ++k) {
      snprintf(format_line, zerork::MAX_FILE_LINE,
               "%7.2f",
               scan_temperature[k]);
      *report += std::string(format_line);

      for(int m=0; m<block_size; m++) {
        list_id = m+block_size*j;
        snprintf(format_line,zerork::MAX_FILE_LINE,
                 " %10.3e",
                 probability->at(list_id).v_double[k]);
        *report += std::string(format_line);
      }
      *report += "\n";
    }

    if(j < num_blocks-1) {
      *report += "\n\n";
    }
  } // for loop-j over the number of complete blocks

  // complete any remaining partial blocks
  // write new block header
  if(rem_columns > 0) {
    *report += "\n\n";
    snprintf(format_line,zerork::MAX_FILE_LINE,
             "# gnuplot data block (index): %d\n",
             num_blocks);
    *report += std::string(format_line);
    *report += "# column 1: [K] scan temperature\n";
    for(int k=0; k<rem_columns; ++k) {
      list_id = k+block_size*num_blocks;
      step_id = probability->at(list_id).id;
      mech->getReactionNameDirOfStep(step_id, &step_name);

      snprintf(format_line,zerork::MAX_FILE_LINE,
               "# column %d: [1/s] step %4d rate {max =%10.3e %s }\n",
               k+2,
               step_id,
               SortableVectorMax(&probability->at(list_id)),
               step_name.c_str());
      *report += std::string(format_line);
    }
    *report += "#------------------------------------------------------------------------------\n";
    // end new block header
    for(int k=0; k<num_scans; ++k) {
      snprintf(format_line, zerork::MAX_FILE_LINE,
               "%7.2f",
               scan_temperature[k]);
      *report += std::string(format_line);

      for(int m=0; m<rem_columns; m++) {
        list_id = m+block_size*num_blocks;
        snprintf(format_line,zerork::MAX_FILE_LINE,
                 " %10.3e",
                 probability->at(list_id).v_double[k]);
        *report += std::string(format_line);
      }
      *report += "\n";
    }
  } // end if(rem_columns > 0)
}
