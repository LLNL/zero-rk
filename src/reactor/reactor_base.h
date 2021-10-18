#ifndef REACTOR_BASE_H_
#define REACTOR_BASE_H_

#include <string>
#include <vector>
#include <map>

#include "reactor_constants.h"

#include "zerork/mechanism.h"

class ReactorBase
{
 public:
  // pure virtual member functions that must be defined in the derivative
  // classes
  virtual ~ReactorBase() {};
  virtual ReactorError GetTimeDerivative(const double reactor_time,
                                         const double state[],
                                         double derivative[]) = 0;

  virtual ReactorError GetJacobian(const double reactor_time,
                                   const double state[],
                                   double jacobian[]) = 0;

  virtual int GetJacobianPattern(int row_id[],
                                 int col_id[]) = 0;

  // common functions to all derived reactor classes
  // basic accessors
  int GetNumStates() const {return num_states_;}
  int GetNumSpecies() const {return num_species_;}
  int GetNumReactions() const {return num_reactions_;}
  int GetNumSteps() const {return num_steps_;}

  const char * GetMechanismName() const
    {return mechanism_name_.c_str();}
  const char * GetThermodynamicsName() const
    {return thermodynamics_name_.c_str();}
  const char * GetParserLogName() const
    {return parser_log_name_.c_str();}

  //  Accessors and mutators for use in the derived functions
  int GetJacobianSize() const {return jacobian_size_;}
  ReactorError SetJacobianSize(const int jacobian_size);

  MatrixType GetMatrixType() const {return matrix_type_;}
  void SetMatrixType(const MatrixType matrix_type)
    {matrix_type_ = matrix_type;}

  const char * GetReactorInfo() const
    {return reactor_info_.c_str();}
  void SetReactorInfo(const std::string &reactor_info)
    {reactor_info_ = reactor_info;}

  // state name utitlies
  int GetIdOfState(const char *state_name) const;
  const char * GetNameOfStateId(const int state_id) const;
  int BuildStateNamesMap(const std::vector<std::string> &state_names);

  // A-Factor sensitivity utilities
  const double * GetAMultipliers() const {return &a_multipliers_.at(0);}
  double GetAMultiplierOfForwardReactionId(const int reaction_id) const;
  ReactorError SetAMultiplierOfForwardReactionId(const int reaction_id,
                                                 const double a_multiplier);
  double GetAMultiplierOfReverseReactionId(const int reaction_id) const;
  ReactorError SetAMultiplierOfReverseReactionId(const int reaction_id,
                                                 const double a_multiplier);
  double GetAMultiplierOfStepId(const int step_id) const;
  ReactorError SetAMultiplierOfStepId(const int step_id,
                                      const double a_multiplier);

  ReactorError BuildMechanism(const char mechanism_name[],
                              const char thermodynamics_name[],
                              const char parser_log_name[]);
  void DestroyMechanism();

  void GetSpeciesHydrogenCount(int num_atoms[]) const
    {mechanism_->getSpeciesHydrogenCount(num_atoms);}
  void GetSpeciesNitrogenCount(int num_atoms[]) const
    {mechanism_->getSpeciesNitrogenCount(num_atoms);}
  void GetSpeciesCarbonCount(int num_atoms[]) const
    {mechanism_->getSpeciesCarbonCount(num_atoms);}
  void GetSpeciesOxygenCount(int num_atoms[]) const
    {mechanism_->getSpeciesOxygenCount(num_atoms);}
  void GetSpeciesMolecularWeight(double molecular_weight[]) const
    {mechanism_->getMolWtSpc(molecular_weight);}

  // TODO: have a constant mechanism that can be returned
  zerork::mechanism * GetMechanism() {return mechanism_;}

  double GetGasConstant() const
    {return mechanism_->getGasConstant();}

 private:

  // data set in BuildMechanism
  int num_species_;
  int num_reactions_;
  int num_steps_;
  std::string mechanism_name_;
  std::string thermodynamics_name_;
  std::string parser_log_name_;
  std::vector<double> a_multipliers_;
  zerork::mechanism *mechanism_;

  // data set in the derived constructor
  int jacobian_size_;
  MatrixType matrix_type_;
  std::string reactor_info_;

  // data set in BuildStateNamesMap
  int num_states_;
  std::map<std::string, int> state_names_map_;
  std::vector<std::string> state_names_;
};



#endif
