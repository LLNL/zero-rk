#ifndef MAGMA_MANAGER_H
#define MAGMA_MANAGER_H

#include <vector>

#include <magma_v2.h>

#include "../hip_la_manager/hip_la_manager.h"

class magma_manager : public hip_la_manager
{
 public:
  magma_manager();
  virtual ~magma_manager();

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

  bool factored_;
  int n_;
  int num_batches_;

  std::vector<double*> data_ptrs_; //CPU data pointers
  std::vector<int> info_;

  //device pointers
  int* info_dev_;
  int* ipiv_dev_;
  int** ipiv_pointers_dev_;
  double* tmp_dev_;
  double* matrix_inverse_dev_;
  double** matrix_inverse_pointers_dev_;
  double** matrix_pointers_dev_;
  double** tmp_pointers_dev_;

  magma_queue_t magma_queue_;
};


#endif
