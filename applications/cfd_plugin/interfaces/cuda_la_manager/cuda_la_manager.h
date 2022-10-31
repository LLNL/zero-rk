#ifndef CUDA_LA_MANAGER_H
#define CUDA_LA_MANAGER_H


class cuda_la_manager
{
 public:
  cuda_la_manager();
  virtual ~cuda_la_manager() {};

  int factor(int num_batches, int n, double* values);
  int solve(int num_batches, int n, const double* rhs, double* soln);

  void set_lu(bool val) {use_lu_ = val;};

 protected:
  bool use_lu_; 

  virtual int factor_invert(int num_batches, int n, double* values) = 0;
  virtual int factor_lu(int num_batches, int n, double* values) = 0;

  virtual int solve_lu(int num_batches, int n, const double* rhs, double* soln) = 0;
  virtual int solve_invert(int num_batches, int n, const double* rhs, double* soln) = 0;
};


#endif
