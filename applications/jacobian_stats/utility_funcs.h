
#ifndef UTILITY_FUNCS_H
#define UTILITY_FUNCS_H

#include <nvector/nvector_serial.h>
#include <string>

void print_sp_matrix(int m,int n, int* colSum, int* rowIdx, double* avals);

double get_val_sp_matrix(int* colSum, int* rowIdx, double* vals, int i, int j);

void permute_sparse_csc(int n,const int* aColSum, const int *aRowIdx, const double *aVals, int* bColSum, int* bRowIdx, double *bVals, int* perm_c);
void permute_sparse_csc_iperm(int n,const int* aColSum, const int *aRowIdx, const double *aVals, int* bColSum, int* bRowIdx, double *bVals, int* iperm_c);

int check_flag(void *flagvalue, const char *funcname, int opt);
void print_state(const int nElem, N_Vector state); 

void insertionSort(int n, double A[]);

const int MAX_NUMBER_STR_LENGTH = 1024;
std::string intToStr(const int i, const char *format);
std::string dblToStr(const double d, const char *format);

const int MAX_REPORT_LINE_LENGTH = 1024;

#endif
