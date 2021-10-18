#include <stdlib.h>
#include <stdio.h>
#include <math.h> // for sqrt if necessary
#include "sparse_matrix_dist.h"


SparseMatrix_dist::SparseMatrix_dist(const int num_rows_loc, const int max_nonzeros_loc, MPI_Comm &comm)
{
  comm_ = comm;
  MPI_Comm_size(comm_, &npes_);
  MPI_Comm_rank(comm_, &my_pe_);

  is_first_factor_      = true;
  num_rows_             = num_rows_loc*npes_;
  num_rows_loc_         = num_rows_loc;
  max_nonzeros_loc_     = max_nonzeros_loc;
  num_nonzeros_loc_     = -1;
  fst_row_              = my_pe_*num_rows_loc_;

  ldb_                  = num_rows_loc;
  ldx_                  = num_rows_loc;

  // Try to have a square processor grid
  double sqroot = sqrt((double) npes_);
  int v = (int) sqroot;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;

  nprow_  = npes_/v;
  npcol_  = v;

  row_sum_            = (int *)malloc(sizeof(int)*(num_rows_loc+1));
  col_id_             = (int *)malloc(sizeof(int)*max_nonzeros_loc);
  rhs_                = (double *)malloc(sizeof(double)*num_rows_loc);
  solution_           = (double *)malloc(sizeof(double)*num_rows_loc);
  matrix_             = (double *)malloc(sizeof(double)*max_nonzeros_loc);

  set_default_options_dist(&options_slu_);
  /* Set the default input options:
        options.Fact = DOFACT;
        options.Equil = YES;
	options.ParSymbFact = NO;
        options.ColPerm = MMD_AT_PLUS_A;
	//if metis is not installed. If so, default is METIS_AT_PLUS_A
        options.RowPerm = LargeDiag;
        options.ReplaceTinyPivot = YES;
        options.Trans = NOTRANS;
        options.IterRefine = SLU_DOUBLE;
        options.SolveInitialized = NO;
        options.RefineInitialized = NO;
	options.num_lookaheads = 10;
	options.lookahead_etree = NO;
	options.SymPattern = NO;
        options.PrintStat = YES;
  */
  options_slu_.IterRefine = NOREFINE;
  options_slu_.lookahead_etree = NO;
  options_slu_.num_lookaheads = 0;
  options_slu_.SymPattern = NO;
  options_slu_.PrintStat = NO;
  options_slu_.ParSymbFact = NO;
  options_slu_.ColPerm = MMD_AT_PLUS_A;
}
SparseMatrix_dist::~SparseMatrix_dist()
{
}

void SparseMatrix_dist::SparseMatrixClean_dist()
{
  if(col_id_ != NULL) {
    SUPERLU_FREE(col_id_);
  }
  if(row_sum_ != NULL) {
    SUPERLU_FREE(row_sum_);
  }
  if(rhs_ != NULL) {
    SUPERLU_FREE(rhs_);
  }
  if(solution_ != NULL) {
    SUPERLU_FREE(solution_);
  }
  if(matrix_ != NULL) {
    SUPERLU_FREE(matrix_);
  }
  if(!is_first_factor_) {
    PStatFree(&stats_slu_);
    dScalePermstructFree(&ScalePermstruct_);
    dDestroy_LU(num_rows_, &grid_, &LUstruct_);
    if (options_slu_.SolveInitialized)
      dSolveFinalize(&options_slu_, &SOLVEstruct_);
    dLUstructFree(&LUstruct_);
    superlu_gridexit(&grid_);
  }
}

void SparseMatrix_dist::SetupFirstFactor_dist()
{
  superlu_gridinit(comm_, nprow_, npcol_, &grid_);

  int iam = grid_.iam;
  if ( !iam ) {
    int v_major, v_minor, v_bugfix;
    superlu_dist_GetVersionNumber(&v_major, &v_minor, &v_bugfix);
    printf("Library version:\t%d.%d.%d\n", v_major, v_minor, v_bugfix);
    printf("Process grid:\t\t%d X %d\n", (int)grid_.nprow, (int)grid_.npcol);
    fflush(stdout);
  }

  dCreate_CompRowLoc_Matrix_dist(&A_,//SuperMatrix
				 num_rows_, //number of global rows
				 num_rows_, //number of global columns
				 max_nonzeros_loc_, //number of local nonzeros
				 num_rows_loc_, //number of local rows
				 fst_row_, //global row number of the first local row
				 matrix_, //vector of local nonzeros values
				 col_id_, //array of column indices of local nonzeros
				 row_sum_, //array of beginning of rows
				 SLU_NR_loc,
				 SLU_D,
				 SLU_GE);

  dScalePermstructInit(num_rows_, num_rows_, &ScalePermstruct_);

  dLUstructInit(num_rows_, &LUstruct_);

  PStatInit(&stats_slu_); // initialize SuperLU statistics

}

int SparseMatrix_dist::FactorSamePattern_dist(const double matrix[])
{
  const int num_nonzeros_loc = num_nonzeros_loc_;
  int error_flag;
  double berr[1];

  if(is_first_factor_) {
    printf("ERROR: SparseMatrix::FactorSamePattern(...)\n");
    printf("       can not be used for the first factorization\n");
    fflush(stdout);
    return -1; // TODO: add error codes
  }

  // You must destroy the L and U before each new factorization, or
  // you will have a memory leak
  dDestroy_LU(num_rows_, &grid_, &LUstruct_);

  options_slu_.Fact = SamePattern;

  nrhs_ = 0; // Factorize, don't solve

  // copy the matrix to factor to the SuperLU internal storage
  for(int j=0; j<num_nonzeros_loc; ++j) {
    matrix_[j] = matrix[j];
  }

  pdgssvx(&options_slu_,
	  &A_,
	  &ScalePermstruct_,
	  rhs_,
	  ldb_,
	  nrhs_,
	  &grid_,
	  &LUstruct_,
	  &SOLVEstruct_,
	  berr,
	  &stats_slu_,
	  &error_flag);

  if(error_flag != 0) {
    return error_flag;
  }
  // flag > 0, singular matrix, zero diagonal at row,col = flag
  // flag < 0, illegal input
  // positive return values

  return 0;
}

int SparseMatrix_dist::FactorNewPatternCCS_dist(const int new_num_nonzeros,
						const int new_col_id[],
						const int new_row_sum[],
						const double matrix[])
{
  const int num_nonzeros_loc = new_num_nonzeros;
  const int row_sum_size = num_rows_loc_ + 1;
  int error_flag;
  double berr[1];

  num_nonzeros_loc_ = num_nonzeros_loc;

  if(is_first_factor_) {
    SetupFirstFactor_dist();
    is_first_factor_ = false;
  } else {
    // You must destroy the L and U before each new factorization, or
    // you will have a memory leak
    dDestroy_LU(num_rows_, &grid_, &LUstruct_);
  }
  options_slu_.Fact = DOFACT;

  nrhs_ = 0;

  // copy the new col_id to the SuperLU internal storage
  for(int j=0; j<num_nonzeros_loc; ++j) {
    col_id_[j] = new_col_id[j];
  }
  // copy the new row_sum to the SuperLU internal storage
  for(int j=0; j<row_sum_size; ++j) {
    row_sum_[j] = new_row_sum[j];
  }

  // update the size of the matrix in the SuperLU structure
  ((NRformat_loc *)A_.Store)->nnz_loc = num_nonzeros_loc;

  // copy the matrix to factor to the SuperLU internal storage
  for(int j=0; j<num_nonzeros_loc; ++j) {
    matrix_[j] = matrix[j];
  }

  pdgssvx(&options_slu_,
	  &A_,
	  &ScalePermstruct_,
	  rhs_,
	  ldb_,
	  nrhs_,
	  &grid_,
	  &LUstruct_,
	  &SOLVEstruct_,
	  berr,
	  &stats_slu_,
	  &error_flag);

  if(error_flag != 0) {
    return error_flag;
  }
  // flag > 0, singular matrix, zero diagonal at row,col = flag
  // flag < 0, illegal input
  // positive return values

  return 0;
}

int SparseMatrix_dist::FactorSamePatternCCS_dist(const int new_num_nonzeros,
						const int new_col_id[],
						const int new_row_sum[],
						const double matrix[])
{
  const int num_nonzeros_loc = new_num_nonzeros;
  const int row_sum_size = num_rows_loc_ + 1;
  int error_flag;
  double berr[1];

  num_nonzeros_loc_ = num_nonzeros_loc;

  if(is_first_factor_) {
    SetupFirstFactor_dist();
    is_first_factor_ = false;
  } else {
    // You must destroy the L and U before each new factorization, or
    // you will have a memory leak
    dDestroy_LU(num_rows_, &grid_, &LUstruct_);
  }
  options_slu_.Fact = SamePattern;

  nrhs_ = 0;

  // copy the new col_id to the SuperLU internal storage
  for(int j=0; j<num_nonzeros_loc; ++j) {
    col_id_[j] = new_col_id[j];
  }
  // copy the new row_sum to the SuperLU internal storage
  for(int j=0; j<row_sum_size; ++j) {
    row_sum_[j] = new_row_sum[j];
  }

  // update the size of the matrix in the SuperLU structure
  ((NRformat_loc *)A_.Store)->nnz_loc = num_nonzeros_loc;

  // copy the matrix to factor to the SuperLU internal storage
  for(int j=0; j<num_nonzeros_loc; ++j) {
    matrix_[j] = matrix[j];
  }

  pdgssvx(&options_slu_,
	  &A_,
	  &ScalePermstruct_,
	  rhs_,
	  ldb_,
	  nrhs_,
	  &grid_,
	  &LUstruct_,
	  &SOLVEstruct_,
	  berr,
	  &stats_slu_,
	  &error_flag);

  if(error_flag != 0) {
    return error_flag;
  }
  // flag > 0, singular matrix, zero diagonal at row,col = flag
  // flag < 0, illegal input
  // positive return values

  return 0;
}

int SparseMatrix_dist::Solve_dist(const double rhs[], double solution[])
{
  const int num_rows_loc = num_rows_loc_;
  int error_flag;
  double berr[1];

  if(is_first_factor_) {
    printf("ERROR: SparseMatrix::Solve(...)\n");
    printf("       can not be used before the first factorization\n");
    fflush(stdout);
    return -1; // TODO: add error codes
  }

  options_slu_.Fact = FACTORED;

  nrhs_ = 1; //Solve

  // copy the right hand side to the SuperLU internal storage
  for(int j=0; j<num_rows_loc; ++j) {
    rhs_[j] = rhs[j];
  }

  pdgssvx(&options_slu_,
	  &A_,
	  &ScalePermstruct_,
	  rhs_,
	  ldb_,
	  nrhs_,
	  &grid_,
	  &LUstruct_,
	  &SOLVEstruct_,
	  berr,
	  &stats_slu_,
	  &error_flag);

  // copy the solution to the output vector
  for(int j=0; j<num_rows_loc; ++j) {
    solution[j] = rhs_[j];
  }
  return error_flag;
}
