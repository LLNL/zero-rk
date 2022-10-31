#include "utility_funcs.h"

#include <stdio.h>
#include <math.h>

#include "zerork/mechanism.h"

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

template <typename T>
std::vector<size_t> sort_indexes_pointer(const size_t n, const T* v) {

  // initialize original index locations
  std::vector<size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

template std::vector<size_t> sort_indexes<int>(const std::vector<int> &v);
template std::vector<size_t> sort_indexes<double>(const std::vector<double> &v);
template std::vector<size_t> sort_indexes_pointer<int>(const size_t n, const int* v);
template std::vector<size_t> sort_indexes_pointer<double>(const size_t n, const double* v);

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

double get_val_sp_matrix_csr(int* rowSum, int* colIdx, double* vals, int i, int j)
{
  int n;
  for(n = rowSum[i]; n < rowSum[i+1]; ++n)
  {
      if(colIdx[n] == j) {return vals[n];}
      else if(colIdx[n] > j) {return 0.0;}
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


void print_sp_matrix_csr(int m, int n, int* rowSum, int* colIdx, double* avals)
{
    int i,j;
    for(i = 0; i<m; ++i)
    {
        for(j = 0; j < n; ++j)
        {
           double val = get_val_sp_matrix_csr(rowSum,colIdx,avals,i,j);
//           printf("%14.8e\t",val); 
           printf("%10.3e\t",val); 
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
  int* perm_c = new int[n];
  for(j = 0; j < n; ++j)
  {
    perm_c[iperm_c[j]] = j;
  }
  permute_sparse_csc(n, aColSum, aRowIdx, aVals, bColSum, bRowIdx, bVals, perm_c);
  delete [] perm_c;
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
int check_cvode_flag(void *flagvalue, const char *funcname, int opt)
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



