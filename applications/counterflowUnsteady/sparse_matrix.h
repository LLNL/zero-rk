#ifndef SPARSE_MATRIX_H_
#define SPARSE_MATRIX_H_

#include <slu_ddefs.h>
#include <slu_util.h>

// Wrapper class to simplify the SuperLU interface

class SparseMatrix
{
 public:
  SparseMatrix(const int num_rows, const int max_nonzeros);
  ~SparseMatrix();

  int FactorNewPatternCCS(const int new_num_nonzeros,
                          const int new_row_id[],
                          const int new_column_sum[],
                          const double matrix[]);

  int FactorSamePattern(const double matrix[]);
  int Solve(const double rhs[], double solution[]);

  bool IsFirstFactor() const {return is_first_factor_;}

  void SetNewNumNonZeros(const int new_num_nonzeros);

 private:
  void SetupFirstFactor();

  bool is_first_factor_;
  int num_rows_;
  int max_nonzeros_;
  int num_nonzeros_;

  SuperMatrix M_slu_;
  SuperMatrix L_slu_;
  SuperMatrix U_slu_;
  SuperMatrix B_slu_;
  SuperMatrix X_slu_;

  char equed_[1];
  int *row_permutation_;
  int *column_permutation_;
  int *column_elim_tree_;    // etree

  double *r_vec_;
  double *c_vec_;

  double *matrix_;
  double *rhs_;
  double *solution_;
  int *row_id_;
  int *column_sum_;

#if SUPERLU_MAJOR_VERSION == 5
  GlobalLU_t Glu_; /* facilitate multiple factorizations with
                         SamePattern_SameRowPerm                  */
#endif
  superlu_options_t options_slu_;
  SuperLUStat_t stats_slu_;
  mem_usage_t  mem_usage_;

};


#endif
