#ifndef HIP_LA_MANAGER_H
#define HIP_LA_MANAGER_H

template<typename T>
class hip_la_manager
{
 public:
  hip_la_manager();
  virtual ~hip_la_manager() {};

  int factor(int num_batches, int n, T* values);
  int solve(int num_batches, int n, const T* rhs, T* soln);

  void set_lu(bool val) {use_lu_ = val;};

 protected:
  bool use_lu_; 

  virtual int factor_invert(int num_batches, int n, T* values) = 0;
  virtual int factor_lu(int num_batches, int n, T* values) = 0;

  virtual int solve_lu(int num_batches, int n, const T* rhs, T* soln) = 0;
  virtual int solve_invert(int num_batches, int n, const T* rhs, T* soln) = 0;
};


#endif
