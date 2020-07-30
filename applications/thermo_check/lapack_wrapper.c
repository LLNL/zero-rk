#include <stdio.h>
#include <stdlib.h>

#include "lapack_wrapper.h"


#ifdef __cplusplus
extern "C" {
#endif

void dgels_(char *trans, int *m, int *n, int *nrhs, double *A, int *lda, 
            double *B, int *ldb, double *work, int *lwork, int *info);

void dgglse_(int *m, int *n, int *p, double *A, int *lda, 
             double *B, int *ldb, double *c, double *d, double *x,
             double *work, int *lwork, int *info);
void dgeev_(char *jobvl, char *jobvr, int *n, double *A, int *lda,
            double *wr, double *wi, double *VL, int *ldvl,
            double *VR, int *ldvr, double *work, int *lwork, int *info);


#ifdef __cplusplus
}
#endif

// wrapper for the LAPACK DGELS function with automatic workspace allocation
// and de-allocation
int LapackDGELS(char transpose,
                int num_rows,
                int num_cols,
                int num_rhs,
                double A[],
                int lda,
                double B[],
                int ldb)
{
  int workspace_size;
  int min_dim;
  int min_size;
  int lapack_info;
  double workspace_query;
  double *workspace = NULL;

  if(transpose == 'n' || transpose == 'N') {
    transpose = 'N';
  } else if(transpose == 't' || transpose == 'T') {
    transpose = 'T';
  } else {
    printf("ERROR: DGELS transpose type %c not recognized\n",transpose);
    fflush(stdout);
    return -1; // consistent with the LAPACK error for an illegal value
               // for the transpose character argument
  } 


  // minimum workspace size
  // see http://www.netlib.org/lapack/explore-html/d8/dde/dgels_8f.html
  // min_dim = minimum matrix dimension
  // min_size = max(1, min_dim + max(min_dim, num_rhs))
  min_dim = ((num_cols < num_rows) ? num_cols : num_rows);
  min_size = min_dim + ((min_dim > num_rhs) ? min_dim : num_rhs);
  min_size = ((min_size > 1) ? min_size : 1);

  workspace_size = -1; // indicates that the workspace size is calculated
                       // and returned in the first element of the 
                       // workspace argument

  // call to determine the workspace size
  dgels_(&transpose,
         &num_rows,
         &num_cols,
         &num_rhs,
         A,
         &lda,
         B,
         &ldb,
         &workspace_query,
         &workspace_size,
         &lapack_info);
  if(lapack_info != 0) {
    printf("ERROR: Workspace size query of DGELS returned failed\n");
    printf("       info value = %d\n",lapack_info);
    fflush(stdout);
    return lapack_info;
  }
  workspace_size = (int)(workspace_query+0.5);
  if(workspace_size < min_size) {
    printf("WARNING: Workspace size query of DGELS returned an optimal\n");
    printf("         size = %d smaller than the minimum dimension %d\n",
           workspace_size, min_size);
    printf("         specified in the LAPACK documentation.\n");
    fflush(stdout);
  }
  workspace = (double *) malloc(sizeof(double)*workspace_size);
  if(workspace == NULL) {
    printf("ERROR: Could not allocate workspace for DGELS with size %d.\n",
           workspace_size);
    fflush(stdout);
    return -9; // consistent with the LAPACK error for an illegal value
               // for the workspace argument
  }
  // call for solution
  dgels_(&transpose,
         &num_rows,
         &num_cols,
         &num_rhs,
         A,
         &lda,
         B,
         &ldb,
         workspace,
         &workspace_size,
         &lapack_info);
  if(lapack_info != 0) {
    printf("ERROR: DGELS solution failed with\n");
    printf("       info value = %d\n",lapack_info);
    fflush(stdout);
   }

  if(workspace != NULL) {
    free(workspace);
  }
 
  return lapack_info;
}

// wrapper for the LAPACK DGGLSE function with automatic workspace allocation
// and de-allocation
int LapackDGGLSE(int num_rows_A,
                 int num_cols_A,
                 int num_rows_B,
                 double A[],
                 int lda,
                 double B[],
                 int ldb,
                 double c[],
                 double d[],
                 double solution_x[])
{
  int workspace_size;
  int min_size;
  int lapack_info;
  double workspace_query;
  double *workspace = NULL;

  // minimum workspace size
  // see http://www.netlib.org/lapack/explore-html/d0/d85/dgglse_8f.html
  // min_size = max(1, num_rows_A+num_cols_A+num_rows_B)
  min_size = num_rows_A+num_cols_A+num_rows_B;
  min_size = ((min_size > 1) ? min_size : 1);

  workspace_size = -1; // indicates that the workspace size is calculated
                       // and returned in the first element of the 
                       // workspace argument

  // call to determine the workspace size
  dgglse_(&num_rows_A,
          &num_cols_A,
          &num_rows_B,
          A,
          &lda,
          B,
          &ldb,
          c,
          d,
          solution_x,
          &workspace_query,
          &workspace_size,
          &lapack_info);

  if(lapack_info != 0) {
    printf("ERROR: Workspace size query of DGGLSE returned failed\n");
    printf("       info value = %d\n",lapack_info);
    fflush(stdout);
    return lapack_info;
  }
  workspace_size = (int)(workspace_query+0.5);
  if(workspace_size < min_size) {
    printf("WARNING: Workspace size query of DGGLSE returned an optimal\n");
    printf("         size = %d smaller than the minimum dimension %d\n",
           workspace_size, min_size);
    printf("         specified in the LAPACK documentation.\n");
    fflush(stdout);
  }
  workspace = (double *) malloc(sizeof(double)*workspace_size);
  if(workspace == NULL) {
    printf("ERROR: Could not allocate workspace for DGGLSE with size %d.\n",
           workspace_size);
    fflush(stdout);
    return -11; // consistent with the LAPACK error for an illegal value
               // for the workspace argument
  }
  // call for solution
  dgglse_(&num_rows_A,
          &num_cols_A,
          &num_rows_B,
          A,
          &lda,
          B,
          &ldb,
          c,
          d,
          solution_x,
          workspace,
          &workspace_size,
          &lapack_info);

  if(lapack_info != 0) {
    printf("ERROR: DGGLSE solution failed with\n");
    printf("       info value = %d\n",lapack_info);
    fflush(stdout);
   }

  if(workspace != NULL) {
    free(workspace);
  }
  return lapack_info;
}

// wrapper for the LAPACK DGEEV function with automatic workspace allocation
// and de-allocation

int LapackDGEEV(char find_left_eigenvectors,
                char find_right_eigenvectors,
                int num_rows,
                double A[],
                int lda,
                double real_eigenvalues[],
                double imag_eigenvalues[],
                double left_eigenvectors[],
                int ldvl,
                double right_eigenvectors[],
                int ldvr)

{
  int workspace_size;
  int min_size;
  int lapack_info;
  double workspace_query;
  double *workspace = NULL;

  if(find_left_eigenvectors == 'n' || find_left_eigenvectors == 'N') {
    find_left_eigenvectors = 'N';
  } else if(find_left_eigenvectors == 'v' || find_left_eigenvectors == 'V') {
    find_left_eigenvectors = 'V';
  } else {
    printf("ERROR: DGEEV find_left_eigenvectors flag %c not recognized.\n",find_left_eigenvectors);
    printf("       Must be 'N' for no left eigenvector computation, or\n");
    printf("       must be 'V' to compute the left eigenvectors.\n");
    fflush(stdout);
    return -1; // consistent with the LAPACK error for an illegal value
               // for the left eigenvector flag
  }
 
  if(find_right_eigenvectors == 'n' || find_right_eigenvectors == 'N') {
    find_right_eigenvectors = 'N';
  } else if(find_right_eigenvectors == 'v' || find_right_eigenvectors == 'V') {
    find_right_eigenvectors = 'V';
  } else {
    printf("ERROR: DGEEV find_right_eigenvectors flag %c not recognized.\n",
           find_right_eigenvectors);
    printf("       Must be 'N' for no right eigenvector computation, or\n");
    printf("       must be 'V' to compute the right eigenvectors.\n");
    fflush(stdout);
    return -2; // consistent with the LAPACK error for an illegal value
               // for the right eigenvector flag
  }
  // minimum workspace size
  // see http://www.netlib.org/lapack/explore-html/d0/d85/dgglse_8f.html
  // min_size = max(1, 3*num_rows) if no eigenvectors are calculated
  //          = max(1, 4*num_rows) if eigenvectors are calculated
  min_size = 3*num_rows;
  if(find_left_eigenvectors == 'V' || find_right_eigenvectors =='V') {
    min_size = 4*num_rows;
  }
  if(min_size < 1) {min_size = 1;}

  workspace_size = -1; // indicates that the workspace size is calculated
                       // and returned in the first element of the 
                       // workspace argument

  // call to determine the workspace size
  dgeev_(&find_left_eigenvectors,
         &find_right_eigenvectors,
         &num_rows,
         A,
         &lda,
         real_eigenvalues,
         imag_eigenvalues,
         left_eigenvectors,
         &ldvl,
         right_eigenvectors,
         &ldvr,
         &workspace_query,
         &workspace_size,
         &lapack_info);

  if(lapack_info != 0) {
    printf("ERROR: Workspace size query of DGEEV returned failed\n");
    printf("       info value = %d\n",lapack_info);
    fflush(stdout);
    return lapack_info;
  }
  workspace_size = (int)(workspace_query+0.5);
  if(workspace_size < min_size) {
    printf("WARNING: Workspace size query of DGEEV returned an optimal\n");
    printf("         size = %d smaller than the minimum dimension %d\n",
           workspace_size, min_size);
    printf("         specified in the LAPACK documentation.\n");
    fflush(stdout);
  }
  workspace = (double *) malloc(sizeof(double)*workspace_size);
  if(workspace == NULL) {
    printf("ERROR: Could not allocate workspace for DGEEV with size %d.\n",
           workspace_size);
    fflush(stdout);
    return -12; // consistent with the LAPACK error for an illegal value
                // for the workspace argument
  }
  // call for solution
  dgeev_(&find_left_eigenvectors,
         &find_right_eigenvectors,
         &num_rows,
         A,
         &lda,
         real_eigenvalues,
         imag_eigenvalues,
         left_eigenvectors,
         &ldvl,
         right_eigenvectors,
         &ldvr,
         workspace,
         &workspace_size,
         &lapack_info);

  if(lapack_info != 0) {
    printf("ERROR: DGEEV solution failed with\n");
    printf("       info value = %d\n",lapack_info);
    fflush(stdout);
   }

  if(workspace != NULL) {
    free(workspace);
  }
  return lapack_info;
}
