#ifndef REACTOR_NVECTOR_SERIAL_CUDA_H
#define REACTOR_NVECTOR_SERIAL_CUDA_H

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "interfaces/hip_la_manager/hip_la_manager.h"
#include "interfaces/hipblas_manager/hipblas_manager.h"
#ifdef ZERORK_HAVE_MAGMA
#include "interfaces/magma_manager/magma_manager.h"
#endif
//#include "interfaces/cusolver_rf_manager/cusolver_rf_manager.h"

#include "zerork/mechanism_cuda.h"

#include "reactor_base.h"

class ReactorNVectorSerialGpu : public ReactorBase
{
 public:
  explicit ReactorNVectorSerialGpu(std::shared_ptr<zerork::mechanism_cuda> mech_ptr);
  virtual ~ReactorNVectorSerialGpu();

  int SetRootTime(double t);
  double GetRootTime();

  N_Vector& GetStateNVectorRef();

  void SetBatchMaskNVector(int reactor_idx, N_Vector batch_mask);

  std::vector<double>& GetReactorWeightsRef() { return weights_; };

#if defined SUNDIALS3 || defined SUNDIALS4
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

  int GetNumStateVariables() { return num_variables_; };

  int GetNumBatchReactors() { return num_reactors_; };
  int GetMinBatchReactors() { return 32; };
  int GetMaxBatchReactors() { return max_num_reactors_; };

  void SetSolveTemperature(bool value);
  void SetStepLimiter(double value);

  void Reset();

 protected:
  void DividedDifferenceJacobian(double t, N_Vector y, N_Vector fy,
                                 double* dense_jacobian);

  int CheckMassFractionsDevice(const double* y);
  int SetTemperatures(const double* scaled_temperatures, double* temperatures);
  int ConcentrationDerivative(const double *inverse_densities, double *ydot);
  int TemperatureDerivative(const double* inverse_densities, const double* y, double* ydot);

  std::shared_ptr<zerork::mechanism_cuda> mech_ptr_;
//  cusolver_rf_manager csrfm_temperature_;
//  cusolver_rf_manager csrfm_no_temperature_;
//  cusolver_rf_manager* csrfm_ptr_;
  std::unique_ptr<hip_la_manager> hip_la_manager_;

  bool solve_temperature_;
  int nnz_;
  int nnz_temperature_;
  int nnz_no_temperature_;

  bool jac_setup_;
  int min_num_reactors_;
  int max_num_reactors_;

  double sqrt_unit_round_;

  std::vector<double> state_data_;
  std::vector<double> tmp1_data_;
  std::vector<double> tmp2_data_;
  std::vector<double> tmp3_data_;

  thrust::device_vector<double> state_data_dev_;
  thrust::device_vector<double> tmp1_data_dev_;
  thrust::device_vector<double> tmp2_data_dev_;
  thrust::device_vector<double> tmp3_data_dev_;
  N_Vector state_;
  N_Vector tmp1_;
  N_Vector tmp2_;
  N_Vector tmp3_;

  double initial_time_;
  double root_time_;

  thrust::host_vector<double> inverse_densities_;
  thrust::host_vector<double> pressures_;
  thrust::host_vector<double> dpdts_;
  thrust::host_vector<double> mean_cx_mass_;

  thrust::host_vector<double> mol_wt_;
  thrust::host_vector<double> inv_mol_wt_;
  thrust::host_vector<double> net_production_rates_;
  thrust::host_vector<double> energy_;
  thrust::host_vector<double> cx_mass_;
  thrust::host_vector<double> forward_rates_of_production_;
  thrust::host_vector<double> creation_rates_;
  thrust::host_vector<double> destruction_rates_;

  thrust::device_vector<double> inverse_densities_dev_;
  thrust::device_vector<double> pressures_dev_;
  thrust::device_vector<double> dpdts_dev_;
  thrust::device_vector<double> initial_temperatures_dev_;
  thrust::device_vector<double> initial_energies_dev_;
  thrust::device_vector<double> mean_cx_mass_dev_;
  thrust::device_vector<double> e_src_dev_;
  thrust::device_vector<double> y_src_dev_;

  thrust::device_vector<double> mol_wt_dev_;
  thrust::device_vector<double> inv_mol_wt_dev_;
  thrust::device_vector<double> net_production_rates_dev_;
  thrust::device_vector<double> energy_dev_;
  thrust::device_vector<double> cx_mass_dev_;
  thrust::device_vector<double> forward_rates_of_production_dev_;
  thrust::device_vector<double> creation_rates_dev_;
  thrust::device_vector<double> destruction_rates_dev_;
  thrust::device_vector<double> concentrations_dev_;
  thrust::device_vector<double> temperatures_dev_;

  thrust::host_vector<double> jacobian_data_;
  thrust::device_vector<double> jacobian_data_dev_;
  thrust::host_vector<int>* jacobian_row_sums_ptr_;
  thrust::host_vector<int> jacobian_row_sums_temperature_;
  thrust::host_vector<int> jacobian_row_sums_no_temperature_;
  thrust::host_vector<int>* jacobian_row_indexes_ptr_;
  thrust::host_vector<int> jacobian_row_indexes_temperature_;
  thrust::host_vector<int> jacobian_row_indexes_no_temperature_;
  thrust::host_vector<int>* jacobian_column_indexes_ptr_;
  thrust::host_vector<int> jacobian_column_indexes_temperature_;
  thrust::host_vector<int> jacobian_column_indexes_no_temperature_;
  thrust::device_vector<int>* jacobian_row_sums_dev_ptr_;
  thrust::device_vector<int> jacobian_row_sums_temperature_dev_;
  thrust::device_vector<int> jacobian_row_sums_no_temperature_dev_;
  thrust::device_vector<int>* jacobian_row_indexes_dev_ptr_;
  thrust::device_vector<int> jacobian_row_indexes_temperature_dev_;
  thrust::device_vector<int> jacobian_row_indexes_no_temperature_dev_;
  thrust::device_vector<int>* jacobian_column_indexes_dev_ptr_;
  thrust::device_vector<int> jacobian_column_indexes_temperature_dev_;
  thrust::device_vector<int> jacobian_column_indexes_no_temperature_dev_;

  thrust::host_vector<double> preconditioner_data_;
  thrust::device_vector<double> preconditioner_data_dev_;
  thrust::host_vector<int> preconditioner_row_sums_;
  thrust::host_vector<int> preconditioner_column_indexes_;
  thrust::device_vector<double>* unit_diagonal_dev_ptr_;
  thrust::device_vector<double> unit_diagonal_temperature_dev_;
  thrust::device_vector<double> unit_diagonal_no_temperature_dev_;
  thrust::device_vector<int> preconditioner_row_sums_dev_;
  thrust::device_vector<int> preconditioner_column_indexes_dev_;

  thrust::host_vector<double> dense_jacobian_;
  thrust::host_vector<double> dense_preconditioner_;
  thrust::device_vector<double> dense_jacobian_dev_;
  thrust::device_vector<double> dense_preconditioner_dev_;

  // non-integer reaction network
  int num_noninteger_jacobian_nonzeros_;
  std::vector<int>* noninteger_sparse_id_ptr_;
  std::vector<int> noninteger_sparse_id_temperature_;
  std::vector<int> noninteger_sparse_id_no_temperature_;
  // Only needed for GPU calculation of non-integer jacobian
  int num_noninteger_jacobian_terms_;
  thrust::device_vector<int>* noninteger_term_id_dev_ptr_;
  thrust::device_vector<int> noninteger_term_id_temperature_dev_;
  thrust::device_vector<int> noninteger_term_id_no_temperature_dev_;
  thrust::device_vector<int> noninteger_concentration_id_dev_;
  thrust::device_vector<int> noninteger_step_id_dev_;
  thrust::device_vector<double> noninteger_multiplier_dev_;
  // Only needed for CPU calculation of non-integer jacobian
  std::vector<double> noninteger_jacobian_;
  thrust::host_vector<double> jacobian_data_nonint_;
  thrust::device_vector<double> jacobian_data_nonint_dev_;

  thrust::device_vector<int> destruction_terms_conc_indexes_dev_;
  thrust::device_vector<int> destruction_terms_reac_indexes_dev_;
  thrust::device_vector<int>* destruction_terms_sparse_indexes_dev_ptr_;
  thrust::device_vector<int> destruction_terms_sparse_indexes_temperature_dev_;
  thrust::device_vector<int> destruction_terms_sparse_indexes_no_temperature_dev_;

  thrust::device_vector<int> creation_terms_conc_indexes_dev_;
  thrust::device_vector<int> creation_terms_reac_indexes_dev_;
  thrust::device_vector<int>* creation_terms_sparse_indexes_dev_ptr_;
  thrust::device_vector<int> creation_terms_sparse_indexes_temperature_dev_;
  thrust::device_vector<int> creation_terms_sparse_indexes_no_temperature_dev_;

  void SetupSparseJacobianArrays();
  int SetupJacobianSparseDevice(realtype t, N_Vector y,N_Vector fy);
  int SparseToDenseDevice(const thrust::device_vector<int>& row_idxs, const thrust::device_vector<int>& col_idxs,
                          const thrust::device_vector<double>& vals, double* dense);
  int FormPreconditionerDevice(double gamma);
  int FactorPreconditioner();
  int AddConstantBlockDiagonal(int n, int nbatch, double constant, double* A_dev);

  std::vector<double> weights_;
  thrust::device_vector<double> step_limiter_;
};

#endif
