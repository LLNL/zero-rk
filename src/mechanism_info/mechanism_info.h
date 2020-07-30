#ifndef MECHANISM_INFO_H_
#define MECHANISM_INFO_H_


#include <string>
#include <vector>

class MechanismInfo
{
 public:
  MechanismInfo(const char mechanism_name[],
                const char thermodynamics_name[],
                const char parser_log_name[]);
  ~MechanismInfo();

  void ConvertMoleToMassFraction(const double mole_fraction[],
                                 double mass_fraction[]) const;
  void ConvertMassToMoleFraction(const double mass_fraction[],
                                 double mole_fraction[]) const;

  int GetNumElements() const;
  int GetNumSpecies() const;
  int GetNumReactions() const;
  int GetNumSteps() const;

  // The dimensions of the reaction stoichiometry matrix are 
  // (# rows, # cols) = (num_species, num_reactions).
  //
  // GetReactionStoichiometryMatrixSize() returns the number of nonzero
  // elements in the reaction stoichiometry matrix.
  int GetReactionStoichiometryMatrixSize() const;
  int GetReactionStoichiometryMatrix(int row_id[], 
                                     int column_id[],
                                     double stoichiometry[]) const;

  void GetNumReactionsPerSpecies(int num_reactions_per_species[]) const;

  // reactor_type = "CV" : only constant volume supported
  // state_type = "YT" : only mass fractions + temperature supported
  // pressure_type = 0 : ignore species influence on pressure dependent
  //                     reactions
  // pressure_type = 1 : only include enhanced species influence on pressure
  //                     dependent reactions
  // pressure_type = 2 : include all species influence on pressure dependent
  //                     reactions
  int GetStateDependencyMatrix(const std::string reactor_type,
                               const std::string state_type,
                               const int pressure_type,
                               std::vector<int> *row_id,
                               std::vector<int> *column_id) const;  

  int GetPressureDependentStep(int is_pressure_dependent_step[]) const;
  int GetPressureDependentReaction(int is_pressure_dependent_reaction[]) const;
  int GetInertSpecies(std::vector<int> *inert_species_id) const;

  const char * GetSpeciesName(const int id) const;
  const char * GetReactionName(const int id) const;

  const char * GetMechanismName() const;
  const char * GetThermodynamicsName() const;
  const char * GetParserLogName() const; 

  const char * GetFileDefinitionOfReaction(const int reaction_id) const;
  const char * GetFileReactionUnits() const;
  const char * GetElementList() const;
  
 private:
  class Impl;
  Impl *impl_;

};

#endif
