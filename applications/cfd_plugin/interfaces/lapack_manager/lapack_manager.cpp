#include <cstring> //std::memcpy
#include <algorithm> //std::max,min

#include "lapack_manager.h"

#ifndef USE_MKL
extern "C"
void dgetrf_(const int*, const int*, double*, const int*, int*, int *);

extern "C"
void dgetrs_(const char* ,const int* ,const int*, double*, const int*, int*, double*, const int*, int*);
#endif

lapack_manager::lapack_manager() :
  factored_(false),
  transpose_(false),
  last_factor_n_(-1)
{
}

lapack_manager::~lapack_manager()
{
}

int lapack_manager::factor(const int m, const int n,
                           const std::vector<double>& matrix)
{
  if(matrix.size() != m*n) {
    return 1;
  }
  return this->factor(m, n, &(matrix[0]));
}

int lapack_manager::factor(const int m, const int n,
                           const double* matrix)
{
  lu_.resize(m*n);
  ipiv_.resize(std::min(m,n));
  for(int i = 0; i < m*n; ++i) {
    lu_[i] = matrix[i];
  }
  int lda = m;
  int flag = 0;
  dgetrf_(&m, &n, &(lu_[0]), &lda, &(ipiv_[0]), &flag);
  factored_=true;
  last_factor_n_ = m;
  return flag;
}


int lapack_manager::solve(const std::vector<double> rhs, std::vector<double>* soln) {
  if(rhs.size() != soln->size()) {
    return 1;
  }
  int n = rhs.size();
  return this->solve(n, &(rhs[0]), &((*soln)[0]));
}

int lapack_manager::solve(int n, const double* rhs, double* soln) {
  if(n != last_factor_n_) {
    return 1;
  }
  //todo transpose
  int lda = n;
  int ldb = n;
  int nrhs = 1;
  char* trans = "N";
  if(transpose_) {
    trans = "T";
  }
  int flag = 0;
  std::memcpy(soln, rhs, sizeof(double)*n);
  dgetrs_(trans, &n, &nrhs, &(lu_[0]), &lda, &(ipiv_[0]), soln, &ldb, &flag);
  return flag;
}

void lapack_manager::SetTranspose(const bool transpose)
{
  transpose_ = transpose;
}


