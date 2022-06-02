#include <string>
#include <vector>
#include <map>

#include "../CKconverter/CKReader.h"
#include <zerork/mechanism.h>
#include "mechanism_info.h"

// implementation class for MechanismInfo
class MechanismInfo::Impl
{
 public:
  Impl(const char mechanism_name[],
       const char thermodynamics_name[],
       const char parser_log_name[]);
  ~Impl();

  const zerork::mechanism * GetMechanism() const {return mechanism_;}
  const char * GetMechanismName() const
  {return mechanism_name_.c_str();}
  const char * GetThermodynamicsName() const
  {return thermodynamics_name_.c_str();}
  const char * GetParserLogName() const
  {return parser_log_name_.c_str();}

  const char * GetFileDefinitionOfReaction(const int reaction_id) const;
  const char * GetFileReactionUnits() const {return reaction_units_.c_str();}
  const char * GetElementList() const {return element_list_.c_str();}
  int GetNumElements() const {return num_elements_;}
  int GetPressureDependentReaction(int is_pressure_dependent_reaction[]) const;

  int GetReactionStoichiometryMatrixSize() const
  {return (int)reaction_stoichiometry_coef_.size();}
  int GetReactionStoichiometryMatrix(int row_id[], 
                                     int column_id[],
                                     double stoichiometry[]) const;
  void GetNumReactionsPerSpecies(int num_reactions_per_species[]) const;
  int GetInertSpecies(std::vector<int> *inert_species_id) const;

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

 private:
  int SetFileDefinitionsFromParser();
  int SetReactionFileDefinition(ckr::CKReader *ckrobj);
  int SetStoichiometryMatrix(ckr::CKReader *ckrobj);
  int SetReactionDependencyMatrix(ckr::CKReader *ckrobj);
  int SetPressureDependentReaction(ckr::CKReader *ckrobj);

  void SetNumReactionsPerSpecies();
  int GetStateDependencyMatrix_CV_YT(const int pressure_type,
                                     std::vector<int> *row_id,
                                     std::vector<int> *column_id) const;

  int num_elements_;
  std::string mechanism_name_;
  std::string thermodynamics_name_;
  std::string parser_log_name_;
  std::string element_list_;
  std::string reaction_units_;
  std::vector<std::string> reaction_file_definition_;

  std::vector<int> reaction_stoichiometry_row_id_;    // species  index
  std::vector<int> reaction_stoichiometry_column_id_; // reaction index
  std::vector<double> reaction_stoichiometry_coef_;
  int num_pressure_dependent_reaction_;
  std::vector<int> is_pressure_dependent_reaction_;
  std::vector<int> num_reactions_per_species_;

  zerork::mechanism *mechanism_;

};

MechanismInfo::Impl::Impl(const char mechanism_name[],
                          const char thermodynamics_name[],
                          const char parser_log_name[])
{
  mechanism_name_      = std::string(mechanism_name);
  thermodynamics_name_ = std::string(thermodynamics_name);
  parser_log_name_     = std::string(parser_log_name);

  reaction_file_definition_.clear();

  mechanism_ = new zerork::mechanism(mechanism_name,
                                     thermodynamics_name,
                                     parser_log_name);

  // reparse to get the mechanism file strings
  SetFileDefinitionsFromParser();
}
MechanismInfo::Impl::~Impl()
{
  delete mechanism_;
}

const char * 
  MechanismInfo::Impl::GetFileDefinitionOfReaction(const int reaction_id) const
{
  if(static_cast<int>(reaction_file_definition_.size()) != 
     mechanism_->getNumReactions()) {

    printf("ERROR: In MechanismInfo::GetFileDefinitionOfReaction(id),\n");
    printf("       file definition list size %d\n",
           static_cast<int>(reaction_file_definition_.size()));
    printf("       does not match the mechanism reaction count %d\n",
           mechanism_->getNumReactions());

    return "reaction file definition list wrong size\n";
  }
  if(reaction_id < 0 || reaction_id >= mechanism_->getNumReactions()) {
    printf("ERROR: In MechanismInfo::GetFileDefinitionOfReaction(id),\n");
    printf("       reaction id = %d is out-of-bounds [0, %d]\n",
           reaction_id, mechanism_->getNumReactions()-1);
    fflush(stdout);
    return "reaction index out-of-bounds\n";
  }

  return reaction_file_definition_[reaction_id].c_str();
}

int MechanismInfo::Impl::SetFileDefinitionsFromParser()
{
  ckr::CKReader *ckrobj;
  ckrobj = new ckr::CKReader();

  if(!(ckrobj->read(mechanism_name_,
                    thermodynamics_name_,
                    ""))) {
    // TODO: add another CKReader object read step to include the parser log
    printf("ERROR: In MechanismInfo::SetFileDefinitionsFromParser(),\n");
    printf("       could not parse mechanism file %s, and\n",
           mechanism_name_.c_str());
    printf("       thermodynamics file %s\n",
           thermodynamics_name_.c_str());
    fflush(stdout);
    return -1;
  }

  // set the pressure dependent reaction list
  SetPressureDependentReaction(ckrobj);

  // set the pressure dependent reaction list
  SetStoichiometryMatrix(ckrobj);

  // set the string vector containing each reaction's file definition
  SetReactionFileDefinition(ckrobj);

  // set the count of reactions in which each species is involved
  SetNumReactionsPerSpecies();

  // set the reaction units string
  reaction_units_.clear();
  switch(ckrobj->units.ActEnergy) {

    case ckr::Cal_per_Mole:
      reaction_units_ = std::string("CAL/MOLE");
      break;
    case ckr::Kcal_per_Mole:
      reaction_units_ = std::string("KCAL/MOLE");
      break;
    case ckr::Joules_per_Mole:
      reaction_units_ = std::string("JOULES/MOLE");
      break;
    case ckr::Kjoules_per_Mole:
      reaction_units_ = std::string("KJOULES/MOLE");
      break;
    case ckr::Kelvin:
      reaction_units_ = std::string("KELVINS");
      break;
    case ckr::Electron_Volts:
      reaction_units_ = std::string("EVOLTS");
      break;
  }
  if(ckrobj->units.Quantity == ckr::Molecules) {
    reaction_units_ += std::string(" MOLECULES");
  }
  
  // set the element list string separated by a single space
  element_list_.clear();
  num_elements_ = static_cast<int>(ckrobj->elements.size());
  for(int j=0; j<num_elements_; ++j) {
    element_list_ += ckrobj->elements[j].name;
    if(j != num_elements_-1) {
      element_list_ += std::string(" ");
    }
  }
  
  delete ckrobj;

  return 0;
}

int
  MechanismInfo::Impl::SetReactionFileDefinition(ckr::CKReader *ckrobj)
{
  std::string reaction_string;

  const int num_reactions = static_cast<int>(ckrobj->reactions.size());

  for(int j=0; j<num_reactions; ++j) {
    const int num_lines = static_cast<int>(ckrobj->reactions[j].lines.size());
    //const int num_comments = 
    //  static_cast<int>(ckrobj->reactions[j].comment.size());

    //printf("! Reaction %d (comments %d)\n",j,num_comments);
    for(int k=0; k<num_lines; ++k) {
      if(k == 0) {
        reaction_string.clear();
      }
      reaction_string += ckrobj->reactions[j].lines[k] + std::string("\n");
      //printf("%s\n",ckrobj->reactions[j].lines[k].c_str());
    }
    reaction_file_definition_.push_back(reaction_string);
    //for(int k=0; k<num_comments; ++k) {
    //  printf("! reaction %d comment %d: %s\n",
    //         j,k,ckrobj->reactions[j].comment[k].c_str());
    //}
  }

  return 0;
}

// The reaction stoichiometry matrix has num_species rows and num_reactions
// columns
int MechanismInfo::Impl::SetStoichiometryMatrix(ckr::CKReader *ckreader)
{
  const int num_reactions = mechanism_->getNumReactions();
  const int num_species   = mechanism_->getNumSpecies();
  int num_nonzero=0;
  std::map<std::string, double> species_coef;
  std::map<std::string, double>::iterator iter;

  reaction_stoichiometry_row_id_.clear();
  reaction_stoichiometry_column_id_.clear();
  reaction_stoichiometry_coef_.clear();

  for(int j=0; j<num_reactions; ++j) {
    
    species_coef.clear();
    // record the stoichiometric coefficients for the reactants
    for(size_t k=0; k<ckreader->reactions[j].reactants.size(); ++k) {

      std::string species_name = ckreader->reactions[j].reactants[k].name;
      double coef = ckreader->reactions[j].reactants[k].number;
      // check if species_name is already in the reaction coefficient list
      iter = species_coef.find(species_name);
      // reactant stoichiometric coefficients are negative
      if(iter == species_coef.end()) {
        // species name not found
        species_coef[species_name] = -coef;  
      } else {
        species_coef[species_name] -= coef;
      }
    }
    // record the stoichiometric coefficients for the products
    for(size_t k=0; k<ckreader->reactions[j].products.size(); ++k) {

      std::string species_name = ckreader->reactions[j].products[k].name;
      double coef = ckreader->reactions[j].products[k].number;
      // check if species_name is already in the reaction coefficient list
      iter = species_coef.find(species_name);
      // product stoichiometric coefficients are positive
      if(iter == species_coef.end()) {
        // species name not found
        species_coef[species_name] = coef;  
      } else {
        species_coef[species_name] += coef;
      }
    }
    // record stoichiometry coefficients for each reaction column
    for(iter = species_coef.begin(); iter != species_coef.end(); ++iter) {
      int species_id = mechanism_->getIdxFromName((iter->first).c_str());
      if(iter->second != 0.0 && species_id >= 0 && species_id < num_species) {

        reaction_stoichiometry_row_id_.push_back(species_id);
        reaction_stoichiometry_column_id_.push_back(j);
        reaction_stoichiometry_coef_.push_back(iter->second);
        ++num_nonzero;
      }
    }
  } // for(int j=0; j<num_reactions; ++j)
  return num_nonzero;
}

int MechanismInfo::Impl::GetReactionStoichiometryMatrix(
    int row_id[], 
    int column_id[],
    double stoichiometry[]) const
{
  const int num_nonzero = GetReactionStoichiometryMatrixSize();
  for(int j=0; j<num_nonzero; ++j) {
    row_id[j]        = reaction_stoichiometry_row_id_[j];
    column_id[j]     = reaction_stoichiometry_column_id_[j];
    stoichiometry[j] = reaction_stoichiometry_coef_[j];
  }
  return num_nonzero;
}

// requires SetReactionStoichiometryMatrix() to be called first
void MechanismInfo::Impl::SetNumReactionsPerSpecies()
{
  const int num_species   = mechanism_->getNumSpecies();
  const int num_nonzero = GetReactionStoichiometryMatrixSize();
 
  num_reactions_per_species_.clear();
  num_reactions_per_species_.assign(num_species, 0);
 
  for(int j=0; j<num_nonzero; ++j) {
    if(reaction_stoichiometry_coef_[j] != 0.0) {
      ++num_reactions_per_species_[reaction_stoichiometry_row_id_[j]];
    }
  }
}

void MechanismInfo::Impl::GetNumReactionsPerSpecies(
    int num_reactions_per_species[]) const
{
  const int num_species   = mechanism_->getNumSpecies();
 
  for(int j=0; j<num_species; ++j) {
    num_reactions_per_species[j] = num_reactions_per_species_[j];
  }
}

int MechanismInfo::Impl::SetPressureDependentReaction(ckr::CKReader *ckreader)
{
  const int num_reactions = mechanism_->getNumReactions();
  is_pressure_dependent_reaction_.assign(num_reactions, 0);
  num_pressure_dependent_reaction_ = 0;

  for(int j=0; j<num_reactions; ++j) {
    
    if(ckreader->reactions[j].kf.type == ckr::PLogInterpolation) {

      is_pressure_dependent_reaction_[j] = 1;
      ++num_pressure_dependent_reaction_;

    } else if(ckreader->reactions[j].isThreeBodyRxn) {

      is_pressure_dependent_reaction_[j] = 1;
      ++num_pressure_dependent_reaction_;

    } else if(ckreader->reactions[j].isFalloffRxn) {
    
      // NOTE the CKReader object may only use "M" to represent a third body
      if(ckreader->reactions[j].thirdBody == std::string("M") ||
         ckreader->reactions[j].thirdBody == std::string("m")) {

        is_pressure_dependent_reaction_[j] = 1;
        ++num_pressure_dependent_reaction_;

      }
    }
  } // for(int j=0; j<num_reactions; ++j)
  return num_pressure_dependent_reaction_;
}
int MechanismInfo::Impl::GetPressureDependentReaction(
    int is_pressure_dependent_reaction[]) const
{
  const int num_reaction = mechanism_->getNumReactions();
  for(int j=0; j<num_reaction; ++j) {
    is_pressure_dependent_reaction[j] = is_pressure_dependent_reaction_[j];
  }
  return num_pressure_dependent_reaction_;
}

int MechanismInfo::Impl::GetInertSpecies(
    std::vector<int> *inert_species_id) const
{
  const int num_species = mechanism_->getNumSpecies();
  inert_species_id->clear();

  for(int j=0; j<num_species; ++j) {

    if(num_reactions_per_species_[j] == 0) {
      inert_species_id->push_back(j);
    }
  }
  return static_cast<int>(inert_species_id->size());
}


int MechanismInfo::Impl::GetStateDependencyMatrix(
    const std::string reactor_type,
    const std::string state_type,
    const int pressure_type,
    std::vector<int> *row_id,
    std::vector<int> *column_id) const
{
  row_id->clear();
  column_id->clear();

  if(reactor_type != "CV" || state_type != "YT") {
    printf("# WARNING: In GetStateDependencyMatrix(...),\n");
    printf("#          reactor_type = %s\n",reactor_type.c_str());
    printf("#          state_type   = %s\n",state_type.c_str());
    printf("#          combination is unsupported.\n");
    printf("#          Supported reactor_type(s): \"CV\"\n");
    printf("#                      state_type(s): \"YT\"\n");
    fflush(stdout);
    return 0;
  }
  return GetStateDependencyMatrix_CV_YT(pressure_type,
                                        row_id,
                                        column_id);
}

int MechanismInfo::Impl::GetStateDependencyMatrix_CV_YT(
    const int pressure_type,
    std::vector<int> *row_id,
    std::vector<int> *column_id) const
{
  const int num_species   = mechanism_->getNumSpecies();
  const int num_reactions = mechanism_->getNumReactions();
  const int num_states    = num_species + 1;
  int cid, rid, global_id;

  std::map<int, int> dependency_matrix;
  std::map<int, int>::iterator dependency_matrix_iter;
  
  std::vector<int> reactant_id, product_id;

  ckr::CKReader *ckreader;
  ckreader = new ckr::CKReader();

  if(!(ckreader->read(mechanism_name_,
                      thermodynamics_name_,
                      ""))) {
    // TODO: add another CKReader object read step to include the parser log
    printf("# ERROR: In MechanismInfo::GetStateDependencyMatrix(),\n");
    printf("#       could not re-parse mechanism file %s, and\n",
           mechanism_name_.c_str());
    printf("#       thermodynamics file %s\n",
           thermodynamics_name_.c_str());
    fflush(stdout);
    return -1;
  }

  dependency_matrix.clear();

  for(int j=0; j<num_reactions; ++j) {

    reactant_id.clear();
    product_id.clear();

    // store the reactants
    for(size_t k=0; k<ckreader->reactions[j].reactants.size(); ++k) {

      std::string species_name = ckreader->reactions[j].reactants[k].name;
      int species_id = mechanism_->getIdxFromName(species_name.c_str());
      if(species_id >= 0 && species_id < num_species) {
        reactant_id.push_back(species_id);
      } else {
        printf("# ERROR: In MechanismInfo::GetStateDependencyMatrix(),\n");
        printf("#        could not find species %s from reaction %d (0-index)\n",species_name.c_str(),j);
        printf("#        file definition:\n");
        printf("%s\n",GetFileDefinitionOfReaction(j));
        fflush(stdout);
        return 0;
      }
    } // for(size_t k=0; k<ckreader->reactions[j].reactants.size(); ++k)

    // store the products
    for(size_t k=0; k<ckreader->reactions[j].products.size(); ++k) {

      std::string species_name = ckreader->reactions[j].products[k].name;
      int species_id = mechanism_->getIdxFromName(species_name.c_str());
      if(species_id >= 0 && species_id < num_species) {
        reactant_id.push_back(species_id);
      } else {
        printf("# ERROR: In MechanismInfo::GetStateDependencyMatrix(),\n");
        printf("#        could not find species %s from reaction %d (0-index)\n",species_name.c_str(),j);
        printf("#        file definition:\n");
        printf("%s\n",GetFileDefinitionOfReaction(j));
        fflush(stdout);
        return 0;
      }
    } // for(size_t k=0; k<ckreader->reactions[j].products.size(); ++k)

    // a perturbation in the reactants influence the net rate for the reactants
    // and the products
    for(size_t k=0; k<reactant_id.size(); ++k) {

      cid = reactant_id[k]; // column index
      // loop over the reactants
      for(size_t m=0; m<reactant_id.size(); ++m) {
        rid = reactant_id[m]; // row index
        global_id = rid + cid*num_states;
        dependency_matrix[global_id] = 1.0;
      }
      // loop over the products
      for(size_t m=0; m<product_id.size(); ++m) {
        rid = product_id[m]; // row index
        global_id = rid + cid*num_states;
        dependency_matrix[global_id] = 1.0;
      }
    } // for(size_t k=0; k<reactant_id.size(); ++k)
 
    // If reversible, then a perturbation in the products influence the 
    // net rate for the reactants and the products
    if(mechanism_->isReversibleReaction(j)) {

      for(size_t k=0; k<product_id.size(); ++k) {
 
        cid = product_id[k]; // column index
        // loop over the reactants
        for(size_t m=0; m<reactant_id.size(); ++m) {
          rid = reactant_id[m]; // row index
          global_id = rid + cid*num_states;
          dependency_matrix[global_id] = 1.0;
        } 
        // loop over the products
        for(size_t m=0; m<product_id.size(); ++m) {
          rid = product_id[m]; // row index
          global_id = rid + cid*num_states;
          dependency_matrix[global_id] = 1.0;
        }
      } // for(size_t k=0; k<product_id.size(); ++k)
    } //if(mechanism_->isReversibleReaction(j))

    // If pressure dependent and including the enhanced 3-body species
    // influence
    if(is_pressure_dependent_reaction_[j] && pressure_type == 1) {

      // a perturbation in the enhanced third bodies will influence
      // the net rate for the reactants and the products
      std::vector<int> enhanced_species_id;
      std::vector<double> alpha; // enhancement factor
      int step_id = mechanism_->getStepIdxOfRxn(j, 1); // forward step
      mechanism_->getEnhancementFactorsOfStep(step_id,
                                              &enhanced_species_id,
                                              &alpha); // alpha is 
                                                       // the enhancement
                                                       // factor minus one
      for(size_t k=0; k<enhanced_species_id.size(); ++k) {
        
        cid = enhanced_species_id[k]; // column index
        if(alpha[k] != 0.0) {
 
          // loop over the reactants
          for(size_t m=0; m<reactant_id.size(); ++m) {
            rid = reactant_id[m]; // row index
            global_id = rid + cid*num_states;
            dependency_matrix[global_id] = 1.0;
          } 
          // loop over the products
          for(size_t m=0; m<product_id.size(); ++m) {
            rid = product_id[m]; // row index
            global_id = rid + cid*num_states;
            dependency_matrix[global_id] = 1.0;
          }
        } // if(alpha[k] != 0.0)

      } // for(size_t k=0; k<enhanced_species_id.size(); ++k)

    } // if(is_pressure_dependent_reaction_[j] && pressure_type == 1)

    // If pressure dependent and including the influence of all species
    // via the pressure or third body concentration in a constant volume
    // reactor
    if(is_pressure_dependent_reaction_[j] && pressure_type == 2) {

      // a perturbation in any species mass fraction will influence
      // the pressure and third-body concentration in a constant volume
      // reactor and thus influence the net rate of the reactants and products
      for(cid=0; cid<num_species; ++cid) {
      
        // loop over the reactants
        for(size_t m=0; m<reactant_id.size(); ++m) {
          rid = reactant_id[m]; // row index
          global_id = rid + cid*num_states;
          dependency_matrix[global_id] = 1.0;
        } 
        // loop over the products
        for(size_t m=0; m<product_id.size(); ++m) {
          rid = product_id[m]; // row index
          global_id = rid + cid*num_states;
          dependency_matrix[global_id] = 1.0;
        }

      } //  for(cid=0; cid<num_species; ++cid)

    } // if(is_pressure_dependent_reaction_[j] && pressure_type == 2)

  } // for(int j=0; j<num_reactions; ++j)
 
  // include the influence of all the species on the temperature as a dense row
  for(int j=0; j<num_states; ++j) {
    rid = num_states-1;
    cid = j;
    global_id = rid + cid*num_states;
    dependency_matrix[global_id] = 1.0;
  }
  // include the influence of the temperature on all the non-inert species 
  // as a dense column
  for(int j=0; j<num_species; ++j) {
    if(num_reactions_per_species_[j] > 0) {
      rid = j;
      cid = num_states-1;
      global_id = rid + cid*num_states;
      dependency_matrix[global_id] = 1.0;
    }
  }

  // iterate through all the map elements, which indicate a dependency
  for(dependency_matrix_iter = dependency_matrix.begin();
      dependency_matrix_iter != dependency_matrix.end(); 
      ++dependency_matrix_iter) {

    global_id = dependency_matrix_iter->first;
    rid = global_id%num_states;
    cid = global_id/num_states;
    row_id->push_back(rid);
    column_id->push_back(cid);
  }

  delete ckreader;
  return static_cast<int>(row_id->size());
}


// ----------------------------------------------------------------------------
// Public facing API
// ----------------------------------------------------------------------------
MechanismInfo::MechanismInfo(const char mechanism_name[],
                             const char thermodynamics_name[],
                             const char parser_log_name[])
{
  impl_ = new Impl(mechanism_name,
                   thermodynamics_name,
                   parser_log_name);
  if(impl_ == NULL) {
    printf("ERROR: MechanismInfo(...) constructor failed\n");
    fflush(stdout);
  } else if(impl_->GetMechanism() == NULL) {
    printf("ERROR: MechanismInfo(...) constructor failed\n");
    fflush(stdout);
  }
}

MechanismInfo::~MechanismInfo()
{
  if(impl_ != NULL) {
    delete impl_;
  }
} 

int MechanismInfo::GetNumElements() const
{
  return impl_->GetNumElements();
}

int MechanismInfo::GetNumSpecies() const
{
  return impl_->GetMechanism()->getNumSpecies();
}

int MechanismInfo::GetNumReactions() const
{
  return impl_->GetMechanism()->getNumReactions();
}

int MechanismInfo::GetNumSteps() const
{
  return impl_->GetMechanism()->getNumSteps();
}

void MechanismInfo::ConvertMoleToMassFraction(const double mole_fraction[],
                                              double mass_fraction[]) const
{
  impl_->GetMechanism()->getYfromX(mole_fraction, mass_fraction);
}
void MechanismInfo::ConvertMassToMoleFraction(const double mass_fraction[],
                                              double mole_fraction[]) const
{
  impl_->GetMechanism()->getXfromY(mass_fraction, mole_fraction);
}

const char * MechanismInfo::GetSpeciesName(const int id) const
{
  return impl_->GetMechanism()->getSpeciesName(id);
}

const char * MechanismInfo::GetReactionName(const int id) const
{
  return impl_->GetMechanism()->getReactionName(id);
}

const char * MechanismInfo::GetMechanismName() const
{
  return impl_->GetMechanismName();
}

const char * MechanismInfo::GetThermodynamicsName() const
{
  return impl_->GetThermodynamicsName();
}

const char * MechanismInfo::GetParserLogName() const
{
  return impl_->GetParserLogName();
}

const char *
  MechanismInfo::GetFileDefinitionOfReaction(const int reaction_id) const
{
  return impl_->GetFileDefinitionOfReaction(reaction_id);
}

const char * MechanismInfo::GetFileReactionUnits() const
{
  return impl_->GetFileReactionUnits();
}
const char * MechanismInfo::GetElementList() const
{
  return impl_->GetElementList();
}

int MechanismInfo::GetPressureDependentReaction(
    int is_pressure_dependent_reaction[]) const
{
  return impl_->GetPressureDependentReaction(is_pressure_dependent_reaction);
}

int MechanismInfo::GetReactionStoichiometryMatrixSize() const
{
  return impl_->GetReactionStoichiometryMatrixSize();
}
int MechanismInfo::GetReactionStoichiometryMatrix(
    int row_id[], 
    int column_id[],
    double stoichiometry[]) const
{
  return impl_->GetReactionStoichiometryMatrix(row_id, 
                                               column_id,
                                               stoichiometry);
}

int MechanismInfo::GetStateDependencyMatrix(
    const std::string reactor_type,
    const std::string state_type,
    const int pressure_type,
    std::vector<int> *row_id,
    std::vector<int> *column_id) const
{
  return impl_->GetStateDependencyMatrix(reactor_type,
                                         state_type,
                                         pressure_type,
                                         row_id,
                                         column_id);
}

void MechanismInfo::GetNumReactionsPerSpecies(
    int num_reactions_per_species[]) const
{
  impl_->GetNumReactionsPerSpecies(num_reactions_per_species);
}

int MechanismInfo::GetInertSpecies(std::vector<int> *inert_species_id) const
{
  return impl_->GetInertSpecies(inert_species_id);
}
