#include <stdlib.h>
#include <stdio.h>

#include "sparse_matrix.h"


SparseMatrix::SparseMatrix(const int num_rows, const int max_nonzeros)
{
  is_first_factor_      = true;
  num_rows_             = num_rows;
  max_nonzeros_         = max_nonzeros;
  num_nonzeros_         = -1;

  row_permutation_    = (int *)malloc(sizeof(int)*num_rows);
  column_permutation_ = (int *)malloc(sizeof(int)*num_rows);
  column_elim_tree_   = (int *)malloc(sizeof(int)*num_rows);
  column_sum_         = (int *)malloc(sizeof(int)*(num_rows+1));
  row_id_             = (int *)malloc(sizeof(int)*max_nonzeros);
  r_vec_              = (double *)malloc(sizeof(double)*num_rows);
  c_vec_              = (double *)malloc(sizeof(double)*num_rows);
  rhs_                = (double *)malloc(sizeof(double)*num_rows);
  solution_           = (double *)malloc(sizeof(double)*num_rows);
  matrix_             = (double *)malloc(sizeof(double)*max_nonzeros);

  set_default_options(&options_slu_);
   /* Set the default input options:
	options.Fact = DOFACT;
        options.Equil = YES;
    	options.ColPerm = COLAMD;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
	options.DiagPivotThresh = 1.0;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES;
   */
  options_slu_.DiagPivotThresh = 0.0;
  options_slu_.ColPerm = MMD_AT_PLUS_A;

}
SparseMatrix::~SparseMatrix()
{
}
void SparseMatrix::SparseMatrixClean()
{
  if(row_permutation_ != NULL) {
    free(row_permutation_);
  }
  if(column_permutation_ != NULL) {
    free(column_permutation_);
  }
  if(row_id_ != NULL) {
    free(row_id_);
  }
  if(column_sum_ != NULL) {
    free(column_sum_);
  }
  if(column_elim_tree_ != NULL) {
    free(column_elim_tree_);
  }
  if(r_vec_ != NULL) {
    free(r_vec_);
  }
  if(c_vec_ != NULL) {
    free(c_vec_);
  }
  if(rhs_ != NULL) {
    free(rhs_);
  }
  if(solution_ != NULL) {
    free(solution_);
  }
  if(matrix_ != NULL) {
    free(matrix_);
  }
  if(!is_first_factor_) {
    StatFree(&stats_slu_);
    Destroy_SuperMatrix_Store(&B_slu_); // SetupFirstFactor
    Destroy_SuperMatrix_Store(&X_slu_); // SetupFirstFactor
    Destroy_SuperMatrix_Store(&M_slu_); // SetupFirstFactor
    Destroy_SuperNode_Matrix(&L_slu_);  // created in dgssvx
    Destroy_CompCol_Matrix(&U_slu_);    // created in dgssvx
  }
}
void SparseMatrix::SetupFirstFactor_row()
{
  dCreate_CompRow_Matrix(&M_slu_,
	                 num_rows_,
			 num_rows_,
			 max_nonzeros_,
			 matrix_,
			 row_id_,
			 column_sum_,
			 SLU_NR,
                         SLU_D,
                         SLU_GE);
  dCreate_Dense_Matrix(&B_slu_,
                       num_rows_,
                       1,
                       rhs_,
                       num_rows_,
                       SLU_DN,
                       SLU_D,
                       SLU_GE);
  dCreate_Dense_Matrix(&X_slu_,
                       num_rows_,
                       1,
                       solution_,
                       num_rows_,
                       SLU_DN,
                       SLU_D,
                       SLU_GE);

   StatInit(&stats_slu_); // initialize SuperLU statistics

}
void SparseMatrix::SetupFirstFactor_col()
{
  dCreate_CompCol_Matrix(&M_slu_,
	                 num_rows_,
			 num_rows_,
			 max_nonzeros_,
			 matrix_,
			 row_id_,
			 column_sum_,
			 SLU_NC,
                         SLU_D,
                         SLU_GE);
  dCreate_Dense_Matrix(&B_slu_,
                       num_rows_,
                       1,
                       rhs_,
                       num_rows_,
                       SLU_DN,
                       SLU_D,
                       SLU_GE);
  dCreate_Dense_Matrix(&X_slu_,
                       num_rows_,
                       1,
                       solution_,
                       num_rows_,
                       SLU_DN,
                       SLU_D,
                       SLU_GE);

   StatInit(&stats_slu_); // initialize SuperLU statistics

}

int SparseMatrix::FactorSamePattern(const double matrix[])
{
  const int num_nonzeros = num_nonzeros_;
  int error_flag;
  int lwork=0;
  void *work=NULL;
  double rpg,rcond;
  double ferr[1],berr[1]; // length is the number of RHS

  if(is_first_factor_) {
//    printf("ERROR: SparseMatrix::FactorSamePattern(...)\n");
//    printf("       can not be used for the first factorization\n");
//    fflush(stdout);
    return -1;
  }

  // You must destroy the L and U before each new factorization, or
  // you will have a memory leak
  Destroy_SuperNode_Matrix(&L_slu_);
  Destroy_CompCol_Matrix(&U_slu_);

  options_slu_.Fact = SamePattern;
  B_slu_.ncol=0; // in dlinsolx1.c example this is supposed to
                 // indicate that a solution is not needed only
                 // the factorization

  // copy the matrix to factor to the SuperLU internal storage
  for(int j=0; j<num_nonzeros; ++j) {
    matrix_[j] = matrix[j];
  }

  dgssvx(&options_slu_,
         &M_slu_,
         column_permutation_,
         row_permutation_,
         column_elim_tree_,
         equed_,
         r_vec_,
         c_vec_,
         &L_slu_,
         &U_slu_,
         work,
         lwork,
         &B_slu_,
         &X_slu_,
         &rpg,
         &rcond,
         ferr,
         berr,
#if SUPERLU_MAJOR_VERSION == 5
         &Glu_,
#endif
         &mem_usage_,
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
int SparseMatrix::FactorNewPatternCCS(const int new_num_nonzeros,
                                      const int new_row_id[],
                                      const int new_column_sum[],
                                      const double matrix[])
{
  const int num_nonzeros = new_num_nonzeros;
  const int column_sum_size = num_rows_ + 1;
  int error_flag;
  int lwork=0;
  void *work=NULL;
  double rpg,rcond;
  double ferr[1],berr[1]; // length is the number of RHS

  num_nonzeros_ = new_num_nonzeros;

  if(is_first_factor_) {

    SetupFirstFactor_col();
    is_first_factor_ = false;

  } else {
    // You must destroy the L and U before each new factorization, or
    // you will have a memory leak
    Destroy_SuperNode_Matrix(&L_slu_);
    Destroy_CompCol_Matrix(&U_slu_);
  }
  options_slu_.Fact = DOFACT;
  B_slu_.ncol=0; // in dlinsolx1.c example this is supposed to
                 // indicate that a solution is not needed only
                 // the factorization

  // copy the new row_id to the SuperLU internal storage
  for(int j=0; j<num_nonzeros; ++j) {
    row_id_[j] = new_row_id[j];
  }
  // copy the new column_sum to the SuperLU internal storage
  for(int j=0; j<column_sum_size; ++j) {
    column_sum_[j] = new_column_sum[j];
  }

  // update the size of the matrix in the SuperLU structure
  ((NCformat *)M_slu_.Store)->nnz = num_nonzeros;

  // copy the matrix to factor to the SuperLU internal storage
  for(int j=0; j<num_nonzeros; ++j) {
    matrix_[j] = matrix[j];
  }

  // Advanced feature that needs non-standard SuperLU options, and additional
  // logic to save the permutations
  //get_perm_c(options_slu_.ColPerm,
  //           &M_slu_,
  //           column_permutation_);

  dgssvx(&options_slu_,
         &M_slu_,
         column_permutation_,
         row_permutation_,
         column_elim_tree_,
         equed_,
         r_vec_,
         c_vec_,
         &L_slu_,
         &U_slu_,
         work,
         lwork,
         &B_slu_,
         &X_slu_,
         &rpg,
         &rcond,
         ferr,
         berr,
#if SUPERLU_MAJOR_VERSION == 5
         &Glu_,
#endif
         &mem_usage_,
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

int SparseMatrix::FactorNewPatternCRS(const int new_num_nonzeros,
                                      const int new_row_id[],
                                      const int new_column_sum[],
                                      const double matrix[])
{
  const int num_nonzeros = new_num_nonzeros;
  const int column_sum_size = num_rows_ + 1;
  int error_flag;
  int lwork=0;
  void *work=NULL;
  double rpg,rcond;
  double ferr[1],berr[1]; // length is the number of RHS

  num_nonzeros_ = new_num_nonzeros;

  if(is_first_factor_) {

    SetupFirstFactor_row();
    is_first_factor_ = false;

  } else {
    // You must destroy the L and U before each new factorization, or
    // you will have a memory leak
    Destroy_SuperNode_Matrix(&L_slu_);
    Destroy_CompCol_Matrix(&U_slu_);
  }
  options_slu_.Fact = DOFACT;
  B_slu_.ncol=0; // in dlinsolx1.c example this is supposed to
                 // indicate that a solution is not needed only
                 // the factorization

  // copy the new row_id to the SuperLU internal storage
  for(int j=0; j<num_nonzeros; ++j) {
    row_id_[j] = new_row_id[j];
  }
  // copy the new column_sum to the SuperLU internal storage
  for(int j=0; j<column_sum_size; ++j) {
    column_sum_[j] = new_column_sum[j];
  }

  // update the size of the matrix in the SuperLU structure
  ((NCformat *)M_slu_.Store)->nnz = num_nonzeros;

  // copy the matrix to factor to the SuperLU internal storage
  for(int j=0; j<num_nonzeros; ++j) {
    matrix_[j] = matrix[j];
  }

  // Advanced feature that needs non-standard SuperLU options, and additional
  // logic to save the permutations
  //get_perm_c(options_slu_.ColPerm,
  //           &M_slu_,
  //           column_permutation_);

  dgssvx(&options_slu_,
         &M_slu_,
         column_permutation_,
         row_permutation_,
         column_elim_tree_,
         equed_,
         r_vec_,
         c_vec_,
         &L_slu_,
         &U_slu_,
         work,
         lwork,
         &B_slu_,
         &X_slu_,
         &rpg,
         &rcond,
         ferr,
         berr,
#if SUPERLU_MAJOR_VERSION == 5
         &Glu_,
#endif
         &mem_usage_,
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

int SparseMatrix::Solve(const double rhs[], double solution[])
{
  const int num_rows = num_rows_;

  int error_flag;
  int lwork=0;
  void *work=NULL;
  double rpg,rcond;
  double ferr[1],berr[1]; // length is the number of RHS

  if(is_first_factor_) {
//    printf("ERROR: SparseMatrix::Solve(...)\n");
//    printf("       can not be used before the first factorization\n");
//    fflush(stdout);
    return -1;
  }

  options_slu_.Fact = FACTORED;
  B_slu_.ncol=1; // in dlinsolx1.c example this is reset to the number of
                 // RHS to indicate that a solution is needed

  // copy the right hand side to the SuperLU internal storage
  for(int j=0; j<num_rows; ++j) {
    rhs_[j] = rhs[j];
  }

  dgssvx(&options_slu_,
         &M_slu_,
         column_permutation_,
         row_permutation_,
         column_elim_tree_,
         equed_,
         r_vec_,
         c_vec_,
         &L_slu_,
         &U_slu_,
         work,
         lwork,
         &B_slu_,
         &X_slu_,
         &rpg,
         &rcond,
         ferr,
         berr,
#if SUPERLU_MAJOR_VERSION == 5
         &Glu_,
#endif
         &mem_usage_,
         &stats_slu_,
         &error_flag);

  // copy the solution to the output vector
  for(int j=0; j<num_rows; ++j) {
    solution[j] = solution_[j];
  }
  return error_flag;
}
