#include <algorithm> //std::max,min
#include <cstring>   //std::memcpy

#include "lapack_manager_z.h"

#ifndef USE_MKL
extern "C" void zgetrf_(const int *, const int *, std::complex<double> *,
                        const int *, int *, int *);

extern "C" void zgetrs_(const char *, const int *, const int *,
                        std::complex<double> *, const int *, int *,
                        std::complex<double> *, const int *, int *);
#endif

lapack_manager_z::lapack_manager_z()
    : factored_(false), transpose_(false), last_factor_n_(-1) {}

lapack_manager_z::~lapack_manager_z() {}

int lapack_manager_z::factor(const int m, const int n,
                             const std::vector<std::complex<double>> &matrix) {
  if (matrix.size() != m * n) {
    return 1;
  }
  return this->factor(m, n, &(matrix[0]));
}

int lapack_manager_z::factor(const int m, const int n,
                             const std::complex<double> *matrix) {
  lu_.resize(m * n);
  ipiv_.resize(std::min(m, n));
  for (int i = 0; i < m * n; ++i) {
    lu_[i] = matrix[i];
  }
  int lda = m;
  int flag = 0;
  zgetrf_(&m, &n, &(lu_[0]), &lda, &(ipiv_[0]), &flag);
  factored_ = true;
  last_factor_n_ = m;
  return flag;
}

int lapack_manager_z::solve(const std::vector<std::complex<double>> rhs,
                            std::vector<std::complex<double>> *soln) {
  if (rhs.size() != soln->size()) {
    return 1;
  }
  int n = rhs.size();
  return this->solve(n, &(rhs[0]), &((*soln)[0]));
}

int lapack_manager_z::solve(int n, const std::complex<double> *rhs,
                            std::complex<double> *soln) {
  if (n != last_factor_n_) {
    return 1;
  }
  int lda = n;
  int ldb = n;
  int nrhs = 1;
  char *trans = "N";
  if (transpose_) {
    trans = "T";
  }
  int flag = 0;
  std::memcpy(soln, rhs, sizeof(std::complex<double>) * n);
  zgetrs_(trans, &n, &nrhs, &(lu_[0]), &lda, &(ipiv_[0]), soln, &ldb, &flag);
  return flag;
}

void lapack_manager_z::SetTranspose(const bool transpose) {
  transpose_ = transpose;
}
