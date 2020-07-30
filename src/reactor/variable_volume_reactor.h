#ifndef VARIABLE_VOLUME_REACTOR_H_
#define VARIABLE_VOLUME_REACTOR_H_

#include "reactor_constants.h" // defines MatrixType, ReactorError

#include "zerork/mechanism.h"

class VariableVolumeReactor
{
 public:
  VariableVolumeReactor(const char mechanism_name[],
                        const char thermodynamics_name[],
                        const char parser_log_name[],
                        const MatrixType matrix_type,
                        Volume *volume_function,
                        HeatLoss *heat_loss_function);
  ~VariableVolumeReactor();
  
  ReactorError GetTimeDerivative(const double reactor_time,
                                 const double state[],
                                 double derivative[]);

  ReactorError GetJacobian(const double reactor_time,
                           const double state[],
                           double jacobian[]);

  ReactorError GetChemicalHeatReleaseRate(const double reactor_time,
                                          const double state[],
                                          double *heat_release_rate);
  int GetJacobianPattern(int row_id[],
                         int col_id[]);

  int GetNumStates() const;
  int GetNumSpecies() const;
  int GetNumReactions() const;
  int GetNumSteps() const;
  int GetJacobianSize() const;
  MatrixType GetMatrixType() const;

  const char * GetMechanismName() const;
  const char * GetThermodynamicsName() const;
  const char * GetParserLogName() const;
  const char * GetReactorInfo() const;
  zerork::mechanism* GetMechanism();

  // species/reaction accessors
  int GetIdOfState(const char *state_name) const;
  const char * GetNameOfStateId(const int state_id) const;

  // A-Factor sensitivity utilities
  double GetAMultiplierOfForwardReactionId(const int reaction_id) const;
  ReactorError SetAMultiplierOfForwardReactionId(const int reaction_id, 
                                                 const double a_multiplier);
  double GetAMultiplierOfReverseReactionId(const int reaction_id) const;
  ReactorError SetAMultiplierOfReverseReactionId(const int reaction_id, 
                                                 const double a_multiplier);
  double GetAMultiplierOfStepId(const int step_id) const;
  ReactorError SetAMultiplierOfStepId(const int step_id, 
                                      const double a_multiplier);

  void SetReferenceScales(const double ref_moles,
                          const double ref_temperature);
  void GetReferenceScales(double *ref_moles,
                          double *ref_temperature) const;

  void SetVolumeMultiplier(const double v_multiplier);
  double GetVolumeMultiplier() const;

  void GetSpeciesHydrogenCount(int num_atoms[]) const;
  void GetSpeciesNitrogenCount(int num_atoms[]) const;
  void GetSpeciesCarbonCount(int num_atoms[]) const;
  void GetSpeciesOxygenCount(int num_atoms[]) const;

  void GetSpeciesMolecularWeight(double molecular_weight[]) const;

  double GetGasConstant() const;

  // specific heat calculation
  double GetMixtureSpecificHeat_Cv(const double state[]);

 private:
  class Impl;
  Impl *impl_;
};


#endif
