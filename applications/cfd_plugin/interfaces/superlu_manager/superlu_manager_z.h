#ifndef SUPERLU_MANAGER_Z_H
#define SUPERLU_MANAGER_Z_H

#include <vector>
#include "slu_zdefs.h"

class superlu_manager_z
{
 public:
  superlu_manager_z();
  virtual ~superlu_manager_z();

  enum matrix_t { CSC, CSR };
  int factor(const std::vector<int>& indexs, const std::vector<int>& sums,
             const std::vector<doublecomplex>& values, superlu_manager_z::matrix_t type);
  int factor(int n, int nnz, const int* indexs, const int* sums,
             const doublecomplex* values, superlu_manager_z::matrix_t type);
  int refactor(const std::vector<doublecomplex>& values);
  int refactor(int nnz, const doublecomplex* values);
  int solve(const std::vector<doublecomplex> rhs, std::vector<doublecomplex>* soln);
  int solve(int n, const doublecomplex* rhs, doublecomplex* soln);

  bool factored() { return this->m_factored; };

 private:
  void setup_memory(int n, int nnz, const int* indexs, const int* sums,
                    const doublecomplex* values);

  SuperMatrix m_M;
  SuperMatrix m_L;
  SuperMatrix m_U;
  SuperMatrix m_B;
  SuperMatrix m_X;
  superlu_options_t m_options;
  SuperLUStat_t m_stats;
  mem_usage_t m_mem_usage;
  char m_equed[1];
  std::vector<double> m_R;
  std::vector<double> m_C;
  std::vector<int> m_rowPermutation;
  std::vector<int> m_colPermutation;
  std::vector<int> m_colElimTree; // etree

  bool m_factored;
  int m_last_factor_n;
  int m_last_factor_nnz;
#if SUPERLU_MAJOR_VERSION == 5
  GlobalLU_t m_Glu;
#endif
};


#endif
