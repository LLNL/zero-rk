
#ifndef UTILITY_FUNCS_H
#define UTILITY_FUNCS_H


double getHighResolutionTime(void);

void print_sp_matrix(int m,int n, int* colSum, int* rowIdx, double* avals);

double get_val_sp_matrix(int* colSum, int* rowIdx, double* vals, int i, int j);

void permute_sparse_csc(int n,const int* aColSum, const int *aRowIdx, const double *aVals, int* bColSum, int* bRowIdx, double *bVals, int* perm_c);
void permute_sparse_csc_iperm(int n,const int* aColSum, const int *aRowIdx, const double *aVals, int* bColSum, int* bRowIdx, double *bVals, int* iperm_c);


#endif
