#ifndef CUBLAS_MANAGER_H
#define CUBLAS_MANAGER_H

#include <vector>

#include "cublas_v2.h"

#include "../cuda_la_manager/cuda_la_manager.h"

class cublas_manager : public cuda_la_manager
{
 public:
  cublas_manager();
  virtual ~cublas_manager();

  bool factored() { return this->factored_; };

 protected:
  int factor_invert(int num_batches, int n, double* values);
  int factor_lu(int num_batches, int n, double* values);
  int solve_invert(int num_batches, int n, const double* rhs, double* soln);
  int solve_lu(int num_batches, int n, const double* rhs, double* soln);

 private:
  void AllocateDeviceMemory();
  void FreeDeviceMemory();
  void setup_memory();
  int cuda_bdmv(int n, int nbatch, double* A_dev, double* B_dev, double* Y_dev);
  void cuda_transpose(double* odata, const double* idata, const int width, const int height);

  cublasHandle_t cublas_handle_;
  cudaError_t cudaStatus_;

  bool factored_;
  int n_;
  int num_batches_;

  std::vector<double*> data_ptrs_; //CPU data pointers
  std::vector<double*> tmp_ptrs_; //CPU data pointers
  std::vector<int> info_;

  //device pointers
  int* info_dev_;
  double* tmp_dev_;
  double* matrix_inverse_dev_;
  double** matrix_inverse_pointers_dev_;
  double** matrix_pointers_dev_;
  double** tmp_pointers_dev_;
};


#endif
