#ifndef COUNTERFLOW_REACTOR_H_
#define COUNTERFLOW_REACTOR_H_

#include "reactor_constants.h" // defines MatrixType, ReactorError

#include "zerork/mechanism.h"

class CounterflowReactor
{
 public:
  CounterflowReactor(const char mechanism_name[],
                     const char thermodynamics_name[],
                     const char parser_log_name[],
                     const MatrixType matrix_type,
                     const double pressure,
                     const bool finite_separation);
  ~CounterflowReactor();

  ReactorError GetTimeDerivative(const double reactor_time,
                                 const double state[],
                                 double derivative[]);

  ReactorError GetTimeDerivativeLimiter(const double reactor_time,
					const double state[],
					const double step_limiter[],
					double derivative[]);

  ReactorError GetTimeDerivativeSteady(const double state[],
                                       const double step_limiter[],
				       double derivative[]);

  ReactorError GetJacobian(const double reactor_time,
                           const double state[],
                           double jacobian[]);

  ReactorError GetJacobianSteady(const double state[],
				 const double convective[],
				 const bool Tfix,
                                 const double ref_momentum,
                                 const double step_limiter[],
				 double jacobian[]);

  ReactorError GetJacobianLimiter(const double reactor_time,
				  const double state[],
				  const double step_limiter[],
				  double jacobian[]);

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

  // State variable name accessors. The state vector has the following names:
  //
  // state[id = 0] = "MassFraction_<first species name id = 0>"
  // state[id = 1] = "MassFraction_<second species name id = 1>"
  // state[id = :]
  // state[id = num_species-1] = "MassFraction_<last species name>"
  // state[id = num_species]   = "RelativeVolume" with units [m^3/kg]
  // state[id = num_species+1] = "Temperature" non-dimensionalized by
  //                              the reference temperature. See
  //                              Set/GetReferenceTemperature().
  // state[id = num_species+2] = "Momentum"
  // state[id = num_species+3] = "PStrain"
  //
  // For example, if molecular hydrogen was the first species defined in the
  // mechanism input file as "H2" then
  //
  // state[id = 0] = "MassFraction_H2"  Note that it matches the case used
  //                                    in the mechanism input file.
  // Note that if the state_name is not found a value of -1 is returned.
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

  void SetReferenceTemperature(const double ref_temperature);
  double GetReferenceTemperature() const;
  void SetPressure(const double pressure);
  double GetPressure() const;

  void GetSpeciesHydrogenCount(int num_atoms[]) const;
  void GetSpeciesNitrogenCount(int num_atoms[]) const;
  void GetSpeciesCarbonCount(int num_atoms[]) const;
  void GetSpeciesOxygenCount(int num_atoms[]) const;

  void GetSpeciesMolecularWeight(double molecular_weight[]) const;

  double GetGasConstant() const;

  // specific heat calculation
  double GetMixtureSpecificHeat_Cp(const double state[]);
  double GetMixtureSpecificHeat_Cp(const double state[],
                                   double *species_cp);
  double GetMixtureSpecificHeat_Cp(const double temperature,
                                   const double state[],
                                   double *species_cp);

 private:
  class Impl;
  Impl *impl_;
};


#endif
