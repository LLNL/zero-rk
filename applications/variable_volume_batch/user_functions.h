#ifndef USER_FUNCTIONS_H_
#define USER_FUNCTIONS_H_

#include <vector>
#include <map>

#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros

#ifdef SUNDIALS2
#include <sundials/sundials_direct.h> // DlsMat
#elif defined SUNDIALS3
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#elif defined SUNDIALS4
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#endif

#include <reactor/variable_volume_reactor.h>

#include "VariableVolumeBatchIFP.h"
#include "sparse_matrix.h"
#include "volume_functions.h"

typedef std::map<std::string, double> CompositionMap;

void NormalizeCompositionMap(CompositionMap *composition);


const double USER_GAS_CONSTANT = 8.3144621e3;  // [J/kmol-k]

int VariableVolumeRHS(realtype t,        // [in] ODE system time
                      N_Vector y,        // [in] ODE state vector
                      N_Vector ydot,     // [out] ODE state derivative
		      void *params);     // [in/out]

int VariableVolumeRHS_Limit1(realtype t,        // [in] ODE system time
                      N_Vector y,        // [in] ODE state vector
                      N_Vector ydot,     // [out] ODE state derivative
		      void *params);     // [in/out]
int VariableVolumeRHS_Tmin(realtype t,        // [in] ODE system time
                      N_Vector y,        // [in] ODE state vector
                      N_Vector ydot,     // [out] ODE state derivative
		      void *params);     // [in/out]


// cvode version specific definition - TODO: remove the need to switch
// function definitions by hand
//int VariableVolumeDenseJacobian(int N,          // [in] ODE system size
//                                realtype t,     // [in] ODE system time
//                                N_Vector y,     // [in] ODE state vector
//                                N_Vector ydot,  // [in] ODE state derivative
//                                DlsMat Jac,     // [out] ODE Jacobian
//		                void *params,   // [in/out]
//                                N_Vector tmp1,  // [out] N-length workspace
//                                N_Vector tmp2,  // [out] N-length workspace
//                                N_Vector tmp3); // [out] N-length workspace

// cvode version specific definition - TODO: remove the need to switch
// function definitions by hand
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
                                N_Vector tmp1,  // [out] N-length workspace
                                N_Vector tmp2,  // [out] N-length workspace
                                N_Vector tmp3); // [out] N-length workspace

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
  explicit UserData(const int task_num,
                    const char input_filename[],
                    const char volume_filename[],
                    const double initial_pressure,
                    const double initial_temperature,
                    const double fuel_fraction,
		    const CompositionMap &oxidizer_map);
  ~UserData();

  VariableVolumeBatchIFP * GetParser() {return parser_;}
  VariableVolumeReactor * GetReactor() {return reactor_;}
  SparseMatrix *GetSparseMatrix() {return sparse_matrix_;}
  int GetTaskNum() {return task_num_;}
  double GetInitialTime() {return initial_time_;}
  void GetInitialState(double state[]) const;

  double GetVolume(const double t,
                   const double state[]) const;
  double GetVolumeRate(const double t,
                       const double state[]) const;
  double GetMinVolume() const {return volume_->GetMinVolume();}
  double GetTimeAtMinVolume() const {return volume_->GetTimeAtMinVolume();}
  double GetTemperature(const double state[]) const;
  double GetTotalMoles(const double state[]) const;
  double GetPressure(const double t,
                     const double state[]) const;
  void GetMoleFractions(const double state[],
                        double mole_fractions[]) const;
  double GetChemicalHeatReleaseRate(const double t,
                                    const double state[]) const;
  void GetTrackedSpeciesIds(std::vector<int> &species_ids) const;
  void GetTrackedMoleFractions(const double state[],
                               std::vector<double> &mole_fractions) const;
  int GetNumStates() const {return reactor_->GetNumStates();}

  int *    GetRowId()         {return &row_id_[0];}
  int *    GetColumnId()      {return &column_id_[0];}
  int *    GetDiagonalId()    {return &diagonal_id_[0];}
  int *    GetColumnSum()     {return &column_sum_[0];}
  double * GetJacobian()      {return &jacobian_[0];}
  double * GetSavedJacobian() {return &saved_jacobian_[0];}

  int num_preconditioner_setups;
  int num_preconditioner_solves;
  int num_new_jacobians;

 private:
  int BuildInitialState();
  VariableVolumeReactor *reactor_;
  VariableVolumeBatchIFP *parser_;
  SparseMatrix *sparse_matrix_;
  VolumeFromFile *volume_;
  HeatLoss *heat_loss_;
  int task_num_;
  double initial_time_;
  double initial_pressure_;
  double initial_temperature_;
  double fuel_fraction_;
  CompositionMap oxidizer_map_;
  CompositionMap fuel_map_;
  std::vector<double> initial_state_;
  std::vector<int> tracked_species_ids_;

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



#endif
