#include <stdio.h>
#include <vector>

#include "utility_funcs.h"
#include <math.h>

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


void print_sp_matrix(FILE* file, int m, int n, int* colSum, int* rowIdx, double* avals)
{
    int i,j;
    for(i = 0; i<m; ++i)
    {
        for(j = 0; j < n; ++j)
        {
           double val = get_val_sp_matrix(colSum,rowIdx,avals,i,j);
           fprintf(file, "%14.8e\t",val); 
        }
        fprintf(file, "\n");
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

