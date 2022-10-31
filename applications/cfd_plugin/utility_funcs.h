
#ifndef UTILITY_FUNCS_H
#define UTILITY_FUNCS_H

#include <vector>
#include <numeric>
#include <algorithm>

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v);

template <typename T>
std::vector<size_t> sort_indexes_pointer(const size_t n, const T* v);

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

//Matrix markit biz
void print_MM_format_CSC(const char* fname, int nSize,
                         int* rowsum, int* colidx, double* avals);

void print_MM_format_CSR(const char* fname, int nSize,
                         int* rowsum, int* colidx, double* avals);


#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

#endif
