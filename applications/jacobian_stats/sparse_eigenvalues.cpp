#include <stdlib.h>
#include <stdio.h>

#include "sparse_eigenvalues.h"


#ifdef __cplusplus
extern "C" {
#endif

void dgeev_(char *jobvl, char *jobvr, int *n, double *A, int *lda,
            double *wr, double *wi, double *VL, int *ldvl,
            double *VR, int *ldvr, double *work, int *lwork, int *info);

#ifdef __cplusplus
}
#endif


int GetEigenvalues(const int num_rows,
                   const int num_nonzeros,
                   const double sparse_matrix[],
                   const int row_id[], 
                   const int column_sum[],
                   double real_eigenvalues[],
                   double imag_eigenvalues[])
{
  // outputs from the lapack call
  int lapack_info;
  double workspace_query;
  double *dense_matrix, *workspace;

  // fixed values
  char find_left_eigenvectors  = 'N';
  char find_right_eigenvectors = 'N';

  int workspace_size = -1; // indicates that the workspace size is calculated
                           // and returned in the first element of the 
                           // workspace argument
  
  // non-constant arguments for the LAPACK call by reference
  double* right_eigenvectors = NULL;
  double* left_eigenvectors = NULL;
  int ld_eigenvectors=num_rows;
  int copy_num_rows = num_rows;

  if(num_rows < 1) {
    printf("# ERROR: In GetEigenvalues(...),\n");
    printf("#        illegal number of rows = %d\n", num_rows);
    fflush(stdout);
    return -3;  // consistent with the LAPACK error for an illegal value
                // for the number of rows argument
  }

  
  dense_matrix = new double[num_rows*num_rows];
  if(dense_matrix == NULL) {
    printf("# ERROR: In GetEigenvalues(...),\n");
    printf("#        failed to allocate internal dense matrix with size\n");
    printf("#        %d x %d\n",num_rows, num_rows);
    fflush(stdout);
    return -4;  // consistent with the LAPACK error for an illegal value
                // for the dense matrix argument
  }
  for(int j=0; j<(num_rows*num_rows); ++j) {
    dense_matrix[j] = 0.0;
  }

  // copy sparse matrix (num_rows == num_columns)
  for(int j=0; j<num_rows; ++j) { // **NOTE** j is the column number

    int num_elements = column_sum[j+1]-column_sum[j];
    int sparse_id = column_sum[j];

    for(int k=0; k<num_elements; ++k) {
      int m = row_id[sparse_id];
      dense_matrix[m+j*num_rows] = sparse_matrix[sparse_id];
      ++sparse_id;
    }
  }

  // perform a workspace size query
  dgeev_(&find_left_eigenvectors,
         &find_right_eigenvectors,
         &copy_num_rows,
         dense_matrix,
         &copy_num_rows,
         real_eigenvalues,
         imag_eigenvalues,
         left_eigenvectors,
         &ld_eigenvectors,
         right_eigenvectors,
         &ld_eigenvectors,
         &workspace_query,
         &workspace_size,
         &lapack_info);

  if(lapack_info != 0) {
    printf("# ERROR: In GetEigenvalues(...),\n");
    printf("#        workspace size query of DGEEV returned failed with\n");
    printf("#        info value = %d\n",lapack_info);
    fflush(stdout);
    delete [] dense_matrix;
    return lapack_info;
  }

  workspace_size = (int)(workspace_query+0.5);
  workspace = new double[workspace_size];
  if(workspace == NULL) {
    printf("# ERROR: In GetEigenvalues(...),\n");
    printf("#        could not allocate workspace for DGEEV with size %d.\n",
            workspace_size);
 
    fflush(stdout);
    return -12; // consistent with the LAPACK error for an illegal value
                // for the workspace argument
  }

  // call for solution
  dgeev_(&find_left_eigenvectors,
         &find_right_eigenvectors,
         &copy_num_rows,
         dense_matrix,
         &copy_num_rows,
         real_eigenvalues,
         imag_eigenvalues,
         left_eigenvectors,
         &ld_eigenvectors,
         right_eigenvectors,
         &ld_eigenvectors,
         workspace,
         &workspace_size,
         &lapack_info);

  if(lapack_info != 0) {
    printf("# ERROR: In GetEigenvalues(...),\n");
    printf("#        DGEEV solution failed with info value = %d\n",
    lapack_info);
    fflush(stdout);
  }

  if(workspace != NULL) {
    delete [] workspace;
  }
  if(dense_matrix != NULL) {
    delete [] dense_matrix;
  }

  return lapack_info; 
}
