#ifndef CUSOLVER_RF_MANAGER_H
#define CUSOLVER_RF_MANAGER_H

#include <vector>

#include "cusolverRf.h"
#include "../superlu_manager/superlu_manager.h"

class cusolver_rf_manager
{
 public:
  cusolver_rf_manager();
  virtual ~cusolver_rf_manager();

  enum matrix_t { CSC, CSR };
  int factor(int num_batches, int n, int nnz, const int* indexs, const int* sums,
             const double* values, cusolver_rf_manager::matrix_t type);
  int refactor(int num_batches, int n, int nnz, const double* values);
  int solve(int num_batches, int n, const double* rhs, double* soln);

  bool factored() { return this->factored_; };

  void reset();

 private:
  void setup_memory();
  void InitializeCusolverRf();
  void AllocateDeviceMemory();
  void FreeDeviceMemory();

  cusolverRfHandle_t cusolverRfHandle_;
  superlu_manager slum_;

  cusolverStatus_t status_;
  cudaError_t cudaStatus_;

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
};


#endif
