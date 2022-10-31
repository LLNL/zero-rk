#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#include "utility_funcs.h"


double get_val_sp_matrix(int* colSum, int* rowIdx, double* vals, int i, int j)
{
  int n;
  for(n = colSum[j]; n < colSum[j+1]; ++n)
  {
      if(rowIdx[n] == i) {return vals[n];}
      else if(rowIdx[n] > i) {return 0.0;}
  }
  return 0.0;
}


void print_sp_matrix(int m, int n, int* colSum, int* rowIdx, double* avals)
{
    int i,j;
    for(i = 0; i<m; ++i)
    {
        for(j = 0; j < n; ++j)
        {
           double val = get_val_sp_matrix(colSum,rowIdx,avals,i,j);
           printf("%14.8e\t",val); 
        }
        printf("\n");
    }
}


void permute_sparse_csc(int n,const int* aColSum, const int *aRowIdx, const double *aVals, int* bColSum, int* bRowIdx, double *bVals, int* perm_c)
{
// * perm_c: Column permutation vector of size n, which defines the 
// *         permutation matrix Pc; perm_c[i] = j means column i of A is 
// *         in position j in A*Pc.
  int j,k,acol,bnnz;
  bnnz = 0;
  bColSum[0] = bnnz;
  for(j = 0; j < n; ++j) //columns of permuted matrix
    {
      bColSum[j+1] = bColSum[j]; //start the counter
      acol = perm_c[j]; //column of original matrix
      for(k = aColSum[acol]; k < aColSum[acol+1]; ++k)
        {
          bRowIdx[bnnz] = aRowIdx[k];
          bVals[bnnz] = aVals[k];
          ++bnnz;
          bColSum[j+1] = bnnz;
        }
    }
}


void permute_sparse_csc_iperm(int n,const int* aColSum, const int *aRowIdx, const double *aVals, int* bColSum, int* bRowIdx, double *bVals, int* iperm_c)
{
  int j;
  std::vector<int> perm_c(n);
  for(j = 0; j < n; ++j)
  {
    perm_c[iperm_c[j]] = j;
  }
  permute_sparse_csc(n, aColSum, aRowIdx, aVals, bColSum, bRowIdx, bVals, &perm_c[0]);
}


/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
void print_state(const int nElem, N_Vector state)
{
  int j;
  for(j=0; j<nElem; j++) {
    printf("y[%4d]: %.18g\n",j,NV_Ith_S(state,j));
  }
  fflush(stdout);
}
// Simple insertion sort, O(n*n) operations.
void insertionSort(int n, double A[])
{
  int j,k;
  int currMinId;
  double currMin;
  for(j=0; j<n-1; j++) {
    currMinId = j;
    currMin   = A[j];
    for(k=j+1; k<n; k++) {

      if(A[k] < currMin) {
        currMinId = k;
        currMin = A[k];
      }
    }
    if(currMinId != j) {
      // swap with position j
      A[currMinId]=A[j];
      A[j]=currMin;
    } // otherwise the current minumum value is in the correct place
  }
}
std::string intToStr(const int i, const char *format)
{
  char cstr[MAX_NUMBER_STR_LENGTH];
  sprintf(cstr,format,i);
  std::string str = cstr;
  return str;
}

std::string dblToStr(const double d, const char *format)
{
  char cstr[MAX_NUMBER_STR_LENGTH];
  sprintf(cstr,format,d);
  std::string str = cstr;
  return str;
}
