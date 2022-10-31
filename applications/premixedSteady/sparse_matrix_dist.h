#ifndef SPARSE_MATRIX_DIST_H_
#define SPARSE_MATRIX_DIST_H_

//#include <slu_ddefs.h>
//#include <slu_util.h>

#include <superlu_dist_config.h>
#include "superlu_ddefs.h"

#ifdef ZERORK_MPI
#include <mpi.h>
#endif
// Wrapper class to simplify the SuperLU interface

class SparseMatrix_dist
{
 public:
  SparseMatrix_dist(const int num_rows_loc, const int max_nonzeros_loc, MPI_Comm &comm);
  ~SparseMatrix_dist();

  MPI_Comm comm_;
  int my_pe_,npes_;

  int FactorNewPatternCCS_dist(const int new_num_nonzeros,
			       const int new_col_id[],
			       const int new_row_sum[],
			       const double matrix[]);

  int FactorSamePatternCCS_dist(const int new_num_nonzeros,
			       const int new_col_id[],
			       const int new_row_sum[],
			       const double matrix[]);

  int FactorSamePattern_dist(const double matrix[]);
  int Solve_dist(const double rhs[], double solution[]);

  bool IsFirstFactor_dist() const {return is_first_factor_;}

  void SparseMatrixClean_dist();

 private:

  void SetupFirstFactor_dist();

  bool is_first_factor_;
  int num_rows_;
  int num_rows_loc_;
  int max_nonzeros_loc_;
  int num_nonzeros_loc_;
  int fst_row_;

  // For Distributed
  SuperMatrix A_;
  dScalePermstruct_t ScalePermstruct_;
  dLUstruct_t LUstruct_;
  dSOLVEstruct_t SOLVEstruct_;
  gridinfo_t grid_;

  int nrhs_, ldb_, ldx_;
  int nprow_, npcol_;

  double *matrix_;
  double *rhs_;
  double *solution_;
  int *row_sum_;
  int *col_id_;

  superlu_dist_options_t options_slu_;
  SuperLUStat_t stats_slu_;

};


#endif
