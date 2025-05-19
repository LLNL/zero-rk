#ifndef REACTOR_NVECTOR_SERIAL_H
#define REACTOR_NVECTOR_SERIAL_H

#include <memory>
#include <complex>
#include "reactor_base.h"
#include "interfaces/lapack_manager/lapack_manager.h"
#include "interfaces/lapack_manager/lapack_manager_z.h"
#include "interfaces/superlu_manager/superlu_manager.h"
#include "interfaces/superlu_manager/superlu_manager_z.h"
#include "zerork/mechanism.h"

class ReactorNVectorSerial : public ReactorBase
{
 public:
  explicit ReactorNVectorSerial(std::shared_ptr<zerork::mechanism> mech_ptr);
  virtual ~ReactorNVectorSerial();

  N_Vector& GetStateNVectorRef();
  void SetBatchMaskNVector(int reactor_idx, N_Vector batch_mask);

  std::vector<double>& GetReactorWeightsRef() { return weights_; };

  void GetAbsoluteToleranceCorrection(N_Vector correction);

#ifdef SUNDIALS2
  int GetJacobianDense(long int N, double t, N_Vector y, N_Vector fy,
                               DlsMat Jac);
#elif defined SUNDIALS3 || defined SUNDIALS4
  int GetJacobianDense(double t, N_Vector y, N_Vector fy,
                               SUNMatrix Jac);
#else
#error "Unsupported SUNDIALS version"
#endif


  int GetJacobianDenseRaw(long int N, double t, N_Vector y, N_Vector fy,
                       double* Jac);

  int JacobianSetup(double t, N_Vector y, N_Vector fy);

  int JacobianFactor(double gamma);

  int JacobianSolve(double t, N_Vector y, N_Vector fy,
                    N_Vector r, N_Vector z);

  int ComplexJacobianFactor(int k, double alpha, double beta);

  int ComplexJacobianSolve(int k, N_Vector ax, N_Vector bx);

  int RootFunction(double t, N_Vector y, double *root_function);

  int SetRootTime(double t);
  double GetRootTime();

  int GetNumStateVariables();

  int GetNumRootFunctions();

  int GetNumBatchReactors() { return 1; };
  int GetMinBatchReactors() { return 1; };
  int GetMaxBatchReactors() { return 1; };

  void SetSolveTemperature(bool value);
  void SetStepLimiter(double value);

  void Reset();

 protected:
  bool solve_temperature_; 

  std::shared_ptr<zerork::mechanism> mech_ptr_;
  N_Vector state_;
  N_Vector tmp1_;
  N_Vector tmp2_;
  N_Vector tmp3_;

  double initial_time_;
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
  std::vector<double> step_limiter_;

 private:
  superlu_manager slum_;
  lapack_manager lpm_;
  std::vector<double> weights_;

  double sqrt_unit_round_;
  double root_time_;

  std::vector<double> state_data_;
  std::vector<double> tmp1_data_;
  std::vector<double> tmp2_data_;
  std::vector<double> tmp3_data_;

  int nnz_;
  int nnz_temperature_;
  int nnz_no_temperature_;
  std::vector<double> jacobian_data_;
  std::vector<int>* jacobian_column_sums_ptr_;
  std::vector<int> jacobian_column_sums_temperature_;
  std::vector<int> jacobian_column_sums_no_temperature_;
  std::vector<int>* jacobian_row_indexes_ptr_;
  std::vector<int> jacobian_row_indexes_temperature_;
  std::vector<int> jacobian_row_indexes_no_temperature_;
  std::vector<int> jacobian_last_row_indexes_;
  std::vector<int> last_row_indexes_;

  std::vector<double> preconditioner_data_;
  std::vector<int> preconditioner_column_sums_;
  std::vector<int> preconditioner_row_indexes_;

  // Complex J for radau
  int ncmplx_jacs_ = (7 - 1)/2;
  std::vector<lapack_manager_z> lpmz_;
  std::vector<superlu_manager_z> slumz_;
  std::vector<std::vector<doublecomplex>> preconditioner_data_z_;
  std::vector<std::vector<int>> preconditioner_column_sums_z_;
  std::vector<std::vector<int>> preconditioner_row_indexes_z_;

  std::vector<double> dense_jacobian_;
  std::vector<double> dense_preconditioner_;

  // non-integer reaction network
  int num_noninteger_jacobian_nonzeros_;
  std::vector<double> noninteger_jacobian_;
  std::vector<int>* noninteger_sparse_id_ptr_;
  std::vector<int> noninteger_sparse_id_temperature_;
  std::vector<int> noninteger_sparse_id_no_temperature_;

  struct reaction_indexes {
    std::vector<int> concentration_indexes;
    std::vector<int> reaction_indexes;
    std::vector<int> sparse_indexes_temperature;
    std::vector<int> sparse_indexes_no_temperature;
  };
  std::vector<int>* destruction_sparse_indexes_ptr_;
  std::vector<int>* creation_sparse_indexes_ptr_;

  reaction_indexes destruction_terms_;
  reaction_indexes creation_terms_;

  void SetupSparseJacobianArrays();
  int SetupJacobianSparse(realtype t, N_Vector y,N_Vector fy);
#ifdef SUNDIALS2
  int SparseToDense(DlsMat Jac);
#elif defined SUNDIALS3 || defined SUNDIALS4
  int SparseToDense(SUNMatrix Jac);
#endif
  int SparseToDense(const std::vector<int>& sums, const std::vector<int>& idxs,
                    const std::vector<double>& vals, std::vector<double>* dense);
  int SparseToDenseComplex(const std::vector<int>& sums, const std::vector<int>& idxs,
                           const std::vector<doublecomplex>& vals,
			   std::vector<std::complex<double>>* dense);
  int FormPreconditioner(double gamma);
  int FactorPreconditioner();

  void DividedDifferenceJacobian(double t, N_Vector y, N_Vector fy,
                                 std::vector<double>* dense_jacobian);

};

#endif
