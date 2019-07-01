
#ifndef UTILITY_FUNCS_H
#define UTILITY_FUNCS_H


double getHighResolutionTime(void);

void print_sp_matrix(int m,int n, int* colSum, int* rowIdx, double* avals);
void print_sp_matrix_csr(int m,int n, int* rowSum, int* colIdx, double* avals);

double get_val_sp_matrix(int* colSum, int* rowIdx, double* vals, int i, int j);
double get_val_sp_matrix_csc(int* rowSum, int* colIdx, double* vals, int i, int j);

void permute_sparse_csc(int n,const int* aColSum, const int *aRowIdx, const double *aVals, int* bColSum, int* bRowIdx, double *bVals, int* perm_c);
void permute_sparse_csc_iperm(int n,const int* aColSum, const int *aRowIdx, const double *aVals, int* bColSum, int* bRowIdx, double *bVals, int* iperm_c);


int check_cvode_flag(void *flagvalue, const char *funcname, int opt);

// Time step estimates
double chemeq2_dt(void *user_data, double T, double pres,
                  double * massFrac);
double cvode_dti(void *user_data, double T, double pres,
                 double dpdt, double dt, double *massFrac);

#endif
