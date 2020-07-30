#ifndef MECHANISM_STATS_H
#define MECHANISM_STATS_H

#include "utilities/sort_vector.h"
#include "zerork/mechanism.h"

const int MAX_LINE_LENGTH=1024;

int getBasicMechanismReport(zerork::mechanism *mech, std::string *report);

int getSpeciesSummaryReport(zerork::mechanism *mech, std::string *report);


int getInertSpeciesIndex(zerork::mechanism *mech,
                         std::vector<int> *species_id);
void getInertSpeciesReport(zerork::mechanism *mech,
                           std::string *report);

int getLoneReactionSpeciesIndex(zerork::mechanism *mech,
                                std::vector<int> *species_id);
void getLoneReactionSpeciesReport(zerork::mechanism *mech,
                                  std::string *report);

void getAdductSpeciesReport(zerork::mechanism *mech,
			    std::string *report);

void getSourceSpeciesReport(zerork::mechanism *mech,
			    std::string *report);
void getSinkSpeciesReport(zerork::mechanism *mech,
			    std::string *report);

void buildSpeciesReactionCount(zerork::mechanism *mech,
                        std::vector<int> *reaction_count);
void buildSpeciesStepCount(zerork::mechanism *mech,
                           std::vector<int> *step_count);
int buildReactionListOfSpecies(const int species_id,
                               zerork::mechanism *mech, 
                               std::vector<int> *reaction_list);

int updateSpeciesActivityOfStep(const int step_id,
                                zerork::mechanism *mech,
                                std::vector<int> *active_species);
int getInactiveSpeciesAndSteps(std::vector<int> *initial_active_species,
                               zerork::mechanism *mech,
                               std::vector<int> *inactive_species,
                               std::vector<int> *inactive_steps);
void getInactiveReport(std::vector<int> *initial_active_species,
                       zerork::mechanism *mech,
                       std::string *report);

void getStepTypeCount(zerork::mechanism *mech,
                      std::vector<int> *step_type_count);

void scanNegativeRate(zerork::mechanism *mech,
                      const double pressure_ref,
                      const double min_temperature,
                      const double max_temperature,
                      const int num_scans,
                      int *min_rate_step_id,
                      double *min_rate_temperature,
                      double *min_rate);

void getNegativeRateReport(zerork::mechanism *mech,
                           const double pressure_ref,
                           const double min_temperature,
                           const double max_temperature,
			   const int num_scans,
                           std::string *report);

double scanMaxUnimolecularRate(zerork::mechanism *mech,
                               const double pressure_ref,
                               const double min_temperature,
                               const double max_temperature,
                               const int num_scans,
                               double *max_rate_temperature,
                               int *max_rate_step_id);
void getSortedUnimolecularRateReport(zerork::mechanism *mech,
                                     const double pressure_ref,
                                     const double min_temperature,
                                     const double max_temperature,
                                     const int num_scans,
                                     const double min_report_rate,
                                     std::string *report);

int buildUnimolecularListOfSteps(zerork::mechanism *mech,
				 std::vector<int> *step_id);
void getMaxUnimolecularRateReport(zerork::mechanism *mech,
                                  const double pressure_ref,
                                  const double min_temperature,
                                  const double max_temperature,
				  const int num_scans,
                                  std::string *report);
void WriteUnimolecularScanInBlocks(zerork::mechanism *mech,
                           const int block_size,
                           const int num_scans,
			   const int num_uni_steps,
                           const double scan_temperature[],
			   const std::vector<zerork::utilities::SortableVector> *probability,
			   std::string *report);
#endif
