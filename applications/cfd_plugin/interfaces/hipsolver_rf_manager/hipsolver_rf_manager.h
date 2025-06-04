#ifndef HIPSOLVER_RF_MANAGER_H
#define HIPSOLVER_RF_MANAGER_H

#include <vector>

#include <hipsolver/hipsolver.h>
#include "../superlu_manager/superlu_manager.h"

class hipsolver_rf_manager
{
 public:
  hipsolver_rf_manager();
  virtual ~hipsolver_rf_manager();

#if HAVE_HIPSOLVERRF
  enum matrix_t { CSC, CSR };
  int factor(int num_batches, int n, int nnz, const int* indexs, const int* sums,
             const double* values, hipsolver_rf_manager::matrix_t type);
  int refactor(int num_batches, int n, int nnz, const double* values);
  int solve(int num_batches, int n, const double* rhs, double* soln);

  bool factored() { return this->factored_; };

  void reset();

 private:
  void setup_memory();
  void InitializeCusolverRf();
  void AllocateDeviceMemory();
  void FreeDeviceMemory();

  hipsolverRfHandle_t hipsolverRfHandle_;
  superlu_manager slum_;

  hipsolverStatus_t status_;
  hipError_t hipStatus_;

  bool factored_;
  int n_;
  int nnz_;
  int num_batches_;

  int* sums_dev_;
  int* indexes_dev_;
  double* data_dev_;
  int* col_permutation_dev_;
  int* row_permutation_dev_;
  double** data_ptrs_dev_;
  double** rhs_ptrs_dev_;
  double* batch_solve_tmp_dev_;
  double* solve_tmp_dev_;
  std::vector<double*> data_ptrs_; //CPU data pointers
#endif
};


#endif
