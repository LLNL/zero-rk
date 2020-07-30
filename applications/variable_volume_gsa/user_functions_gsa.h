#ifndef USER_FUNCTIONS_H_
#define USER_FUNCTIONS_H_

#include <vector>

#include <cvodes/cvodes.h>            // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros

#ifdef SUNDIALS2
#include <sundials/sundials_direct.h> // DlsMat
#else
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#endif

#include <reactor/variable_volume_reactor.h>

#include "VariableVolumeGSAIFP.h"
#include "sparse_matrix.h"
#include "volume_functions.h"

const double USER_GAS_CONSTANT = 8.3144621e3;  // [J/kmol-k]

int VariableVolumeRHS(realtype t,        // [in] ODE system time
                      N_Vector y,        // [in] ODE state vector
                      N_Vector ydot,     // [out] ODE state derivative
		      void *params);     // [in/out]
#if defined SUNDIALS2
int VariableVolumeDenseJacobian(long int N,     // [in] ODE system size
                                realtype t,     // [in] ODE system time
                                N_Vector y,     // [in] ODE state vector
                                N_Vector ydot,  // [in] ODE state derivative
                                DlsMat Jac,     // [out] ODE Jacobian
		                void *params,   // [in/out]
                                N_Vector tmp1,  // [out] N-length workspace
                                N_Vector tmp2,  // [out] N-length workspace
                                N_Vector tmp3); // [out] N-length workspace

int VariableVolumePreconditionerSetup(realtype t,// [in] ODE system time
                                N_Vector y,      // [in] ODE state vector
                                N_Vector ydot,   // [in] ODE state derivative
                                booleantype jok,
				booleantype *new_j,
				realtype gamma,
		                void *params,    // [in/out]
                                N_Vector tmp1,   // [out] N-length workspace
                                N_Vector tmp2,   // [out] N-length workspace
                                N_Vector tmp3);  // [out] N-length workspace

int VariableVolumePreconditionerSolve(realtype t,// [in] ODE system time
                                N_Vector y,      // [in] ODE state vector
                                N_Vector ydot,   // [in] ODE state derivative
                                N_Vector r,      // [in] jacobian rhs
                                N_Vector z,      // [out]
				realtype gamma,
				realtype delta,
				int lr,
		                void *params,    // [in/out]
				N_Vector tmp);   // [out] N-length workspace
#elif defined SUNDIALS3 || defined SUNDIALS4
int VariableVolumeDenseJacobian(realtype t,     // [in] ODE system time
                                N_Vector y,     // [in] ODE state vector
                                N_Vector ydot,  // [in] ODE state derivative
                                SUNMatrix Jac,     // [out] ODE Jacobian
		                void *params,   // [in/out]
                                N_Vector tmp1,   // [out] N-length workspace
                                N_Vector tmp2,   // [out] N-length workspace
                                N_Vector tmp3);  // [out] N-length workspace

int VariableVolumePreconditionerSetup(realtype t,// [in] ODE system time
                                      N_Vector y,      // [in] ODE state vector
                                      N_Vector ydot,   // [in] ODE state derivative
                                      booleantype jok,
                                      booleantype *new_j,
                                      realtype gamma,
                                      void *params);    // [in/out]

int VariableVolumePreconditionerSolve(realtype t,// [in] ODE system time
                                N_Vector y,      // [in] ODE state vector
                                N_Vector ydot,   // [in] ODE state derivative
                                N_Vector r,      // [in] jacobian rhs
                                N_Vector z,      // [out]
				realtype gamma,
				realtype delta,
				int lr,
                                void *params);    // [in/out]
#endif


class UserData
{
 public:
  explicit UserData(const char filename[], const bool write_parser_log);
  ~UserData();

  VariableVolumeGSAIFP * GetParser() {return parser_;}
  VariableVolumeReactor * GetReactor() {return reactor_;}
  SparseMatrix *GetSparseMatrix() {return sparse_matrix_;}
  double GetInitialTime() {return initial_time_;}
  void GetInitialState(double state[]) const;

  double GetVolume(const double t,
                   const double state[]) const;
  double GetVolumeRate(const double t,
                       const double state[]) const;
  double GetTemperature(const double state[]) const;
  double GetTotalMoles(const double state[]) const;
  double GetPressure(const double t,
                     const double state[]) const;
  void GetMoleFractions(const double state[],
                        double mole_fractions[]) const;
  void GetTrackedSpeciesIds(std::vector<int> &species_ids) const;
  void GetTrackedMoleFractions(const double state[],
                               std::vector<double> &mole_fractions) const;
  int GetNumStates() const {return reactor_->GetNumStates();}

  void SetReferenceHeatRelease(double x) {reference_heat_release_ = x;}
  double GetReferenceHeatRelease() const {return reference_heat_release_;}

  int *    GetRowId()         {return &row_id_[0];}
  int *    GetColumnId()      {return &column_id_[0];}
  int *    GetDiagonalId()    {return &diagonal_id_[0];}
  int *    GetColumnSum()     {return &column_sum_[0];}
  double * GetJacobian()      {return &jacobian_[0];}
  double * GetSavedJacobian() {return &saved_jacobian_[0];}

  int num_preconditioner_setups;
  int num_preconditioner_solves;
  int num_new_jacobians;

  // heat release rate quadrature variables
  std::vector<double> heat_release_rate_powers_;
  int num_quadratures_;
  N_Vector quadrature_;
  double quadrature_rtol_;
  double quadrature_atol_;
  bool quadrature_controls_step_;


 private:
  int BuildInitialState();
  VariableVolumeReactor *reactor_;
  VariableVolumeGSAIFP *parser_;
  SparseMatrix *sparse_matrix_;
  Volume *volume_;
  HeatLoss *heat_loss_;
  double initial_time_;
  std::vector<double> initial_state_;
  std::vector<int> tracked_species_ids_;

  double reference_heat_release_;

  double reference_moles_;

  // TODO: find a better spot for the preconditioner workspace
  //       arrays
  std::vector<int> row_id_;
  std::vector<int> column_id_;
  std::vector<int> diagonal_id_;
  std::vector<int> column_sum_;
  std::vector<double> jacobian_;
  std::vector<double> saved_jacobian_;


};

int ChemicalHRRIntegrand(realtype t,
                         N_Vector y,
                         N_Vector qdot,
                         void *params);

double ChemicalHeatReleaseRate(const double t, const double y[], void *params);
double ChemicalHeatReleaseRatePow2(const double t, const double y[], void *params);
double ChemicalHeatReleaseRatePow4(const double t, const double y[], void *params);
double ChemicalHeatReleaseRatePow6(const double t, const double y[], void *params);
double ChemicalHeatReleaseRatePow8(const double t, const double y[], void *params);
double ChemicalHeatReleaseRatePow10(const double t, const double y[], void *params);

#endif
