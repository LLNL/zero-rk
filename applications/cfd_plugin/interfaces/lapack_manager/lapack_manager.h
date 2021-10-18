#ifndef LAPACK_MANAGER_H
#define LAPACK_MANAGER_H

#include <vector>

class lapack_manager
{
 public:
  lapack_manager();
  virtual ~lapack_manager();

  int factor(int m, int n, const std::vector<double>& matrix);
  int factor(int m, int n, const double* matrix);
  int solve(const std::vector<double> rhs, std::vector<double>* soln);
  int solve(int n, const double* rhs, double* soln);

  bool factored() { return this->factored_; };

  void SetTranspose(const bool transpose);
  void SetLayout(const int layout);

 private:
  std::vector<int> ipiv_;
  std::vector<double> lu_;
  bool factored_;
  bool transpose_;
  int last_factor_n_;
};


#endif
