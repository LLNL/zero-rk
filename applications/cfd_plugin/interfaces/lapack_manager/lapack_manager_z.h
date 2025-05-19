#ifndef LAPACK_MANAGER_Z_H
#define LAPACK_MANAGER_Z_H

#include <complex>
#include <vector>

class lapack_manager_z {
public:
  lapack_manager_z();
  virtual ~lapack_manager_z();

  int factor(int m, int n, const std::vector<std::complex<double>> &matrix);
  int factor(int m, int n, const std::complex<double> *matrix);
  int solve(const std::vector<std::complex<double>> rhs,
            std::vector<std::complex<double>> *soln);
  int solve(int n, const std::complex<double> *rhs, std::complex<double> *soln);

  bool factored() { return this->factored_; };

  void SetTranspose(const bool transpose);
  void SetLayout(const int layout);

private:
  std::vector<int> ipiv_;
  std::vector<std::complex<double>> lu_;
  bool factored_;
  bool transpose_;
  int last_factor_n_;
};

#endif
