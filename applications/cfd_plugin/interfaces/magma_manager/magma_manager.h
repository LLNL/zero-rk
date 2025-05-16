#ifndef MAGMA_MANAGER_H
#define MAGMA_MANAGER_H

#include <vector>

#include <magma_v2.h>

#include "../cuda_la_manager/cuda_la_manager.h"

template<typename T>
class magma_manager : public cuda_la_manager<T>
{
 public:
  magma_manager();
  virtual ~magma_manager();

  bool factored() { return this->factored_; };

 protected:
  int factor_invert(int num_batches, int n, T* values);
  int factor_lu(int num_batches, int n, T* values);
  int solve_invert(int num_batches, int n, const T* rhs, T* soln);
  int solve_lu(int num_batches, int n, const T* rhs, T* soln);

 private:
  void AllocateDeviceMemory();
  void FreeDeviceMemory();
  void setup_memory();
  int cuda_bdmv(int n, int nbatch, T* A_dev, T* B_dev, T* Y_dev);
  void cuda_transpose(T* odata, const T* idata, const int width, const int height);

  void getrf_batched();
  void getri_batched();
  void getrs_batched();

  bool factored_;
  int n_;
  int num_batches_;

  std::vector<T*> data_ptrs_; //CPU data pointers
  std::vector<int> info_;

  //device pointers
  int* info_dev_;
  int* ipiv_dev_;
  int** ipiv_pointers_dev_;
  T* tmp_dev_;
  T* matrix_inverse_dev_;
  T** matrix_inverse_pointers_dev_;
  T** matrix_pointers_dev_;
  T** tmp_pointers_dev_;

  magma_queue_t magma_queue_;
};


#endif
