#ifndef REACTOR_NVECTOR_SERIAL_H
#define REACTOR_NVECTOR_SERIAL_H

#include <memory>
#include "reactor_base.h"
#include "interfaces/lapack_manager/lapack_manager.h"
#include "interfaces/superlu_manager/superlu_manager.h"
#include "zerork/mechanism.h"

class ReactorNVectorSerial : public ReactorBase
{
 public:
  ReactorNVectorSerial(std::shared_ptr<zerork::mechanism> mech_ptr);
  virtual ~ReactorNVectorSerial();

  N_Vector& GetStateNVectorRef();
  void SetBatchMaskNVector(int reactor_idx, N_Vector batch_mask);

  std::vector<double>& GetReactorWeightsRef() { return weights_; };

  void GetAbsoluteToleranceCorrection(N_Vector correction);

#ifdef SUNDIALS2
  int GetJacobianDense(long int N, double t, N_Vector y, N_Vector fy,
                               DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#elif defined SUNDIALS3 || defined SUNDIALS4
  int GetJacobianDense(double t, N_Vector y, N_Vector fy,
                               SUNMatrix Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#else
#error "Unsupported SUNDIALS version"
#endif

  int JacobianSetup(double t, N_Vector y, N_Vector fy,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  int JacobianFactor(double gamma);

  int JacobianSolve(double t, N_Vector y, N_Vector fy,
                    N_Vector r, N_Vector z, N_Vector tmp);

  int RootFunction(double t, N_Vector y, double *root_function);

  int GetNumStateVariables();

  int GetNumRootFunctions();

  int GetNumBatchReactors() { return 1; };
  int GetMinBatchReactors() { return 1; };
  int GetMaxBatchReactors() { return 1; };

  void Reset();

 protected:
  std::shared_ptr<zerork::mechanism> mech_ptr_;
  N_Vector state_;
  N_Vector batch_mask_;

  double pressure_;
  double inverse_density_;
  double dpdt_;
  double e_src_;
  const double* y_src_;
  double initial_temperature_;
  double initial_energy_;
  double mean_cx_mass_;

  std::vector<double> mol_wt_;
  std::vector<double> inv_mol_wt_;
  std::vector<double> net_production_rates_;
  std::vector<double> energy_;
  std::vector<double> cx_mass_;
  std::vector<double> forward_rates_of_production_;
  std::vector<double> creation_rates_;
  std::vector<double> destruction_rates_;
  std::vector<double> concentrations_;

 private:
  superlu_manager slum_;
  lapack_manager lpm_;
  std::vector<double> weights_;

  double sqrt_unit_round_;

  int nnz_;
  std::vector<double> jacobian_data_;
  std::vector<int> jacobian_column_sums_;
  std::vector<int> jacobian_row_indexes_;
  std::vector<int> jacobian_diagonal_indexes_;
  std::vector<int> jacobian_last_row_indexes_;
  std::vector<int> diagonal_indexes_;
  std::vector<int> last_row_indexes_;

  std::vector<double> preconditioner_data_;
  std::vector<int> preconditioner_column_sums_;
  std::vector<int> preconditioner_row_indexes_;

  std::vector<double> dense_jacobian_;
  std::vector<double> dense_preconditioner_;

  // non-integer reaction network
  int num_noninteger_jacobian_nonzeros_;
  std::vector<double> noninteger_jacobian_;
  std::vector<int> noninteger_sparse_id_;

  struct reaction_indexes {
    int concentration_index;
    int reaction_index;
    int sparse_index;
  };

  std::vector<reaction_indexes> destruction_terms_;
  std::vector<reaction_indexes> creation_terms_;

  void SetupSparseJacobianArrays();
  int SetupJacobianSparse(realtype t, N_Vector y,N_Vector fy,
                          N_Vector tmp1,N_Vector tmp2,N_Vector tmp3);
#ifdef SUNDIALS2
  int SparseToDense(DlsMat Jac);
#elif defined SUNDIALS3 || defined SUNDIALS4
  int SparseToDense(SUNMatrix Jac);
#endif
  int SparseToDense(const std::vector<int> sums, const std::vector<int> idxs,
                    const std::vector<double> vals, std::vector<double>* dense);
  int FormPreconditioner(double gamma);
  int FactorPreconditioner();

  void DividedDifferenceJacobian(double t, N_Vector y, N_Vector fy,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3,
                                 std::vector<double>* dense_jacobian);

};

#endif
