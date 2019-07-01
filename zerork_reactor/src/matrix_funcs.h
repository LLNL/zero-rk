#ifndef MATRIX_FUNCS_H
#define MATRIX_FUNCS_H

#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros
#include <sundials/sundials_direct.h> // dense DlsMat types, fcts., and macros

int jac_full_prec_setup(realtype t, N_Vector y, N_Vector fy,
			booleantype jok, booleantype *jcurPtr,
			realtype gamma, void *user_data,
			N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int jac_full_prec_solveV3(realtype t, N_Vector y, N_Vector fy,
			  N_Vector r, N_Vector z, realtype gamma,
			  realtype delta, int lr, void *user_data,
			  N_Vector tmp);

void setupJacobianSparse(realtype t, N_Vector y, N_Vector fy,
                         void* user_data, N_Vector tmp1,
                         N_Vector tmp2, N_Vector tmp3);


void sparseOffDiagThreshCopyGamma(const double tol, const int nsize, const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int *nnzB, double B[], 
			     int BrowId[], int BcolSum[], const double gamma);

void sparseOffDiagThreshCopyGamma_rnorm(const double tol, const int nsize, 
			     const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int *nnzB, double B[], 
				   int BrowId[], int BcolSum[], const double gamma);

void sparseOffDiagThreshCopyGamma_cnorm(const double tol, const int nsize, 
			     const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int *nnzB, double B[], 
				   int BrowId[], int BcolSum[], const double gamma);

void sparseOffDiagThreshCopyGamma_rcmin(const double tol, const int nsize, 
			     const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int *nnzB, double B[], 
				   int BrowId[], int BcolSum[], const double gamma);

void sparseOffDiagMagCopyGamma(const double frac, const int nsize, const int nnzA,
			  const double A[], const int ArowId[],
			  const int AcolSum[], int *nnzB, double B[], 
			  int BrowId[], int BcolSum[], const double gamma);

void sparseUpdateGamma(const int nsize, const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int nnzB, double B[], 
			     int BrowId[], int BcolSum[], const double gamma);

int checkPattern(const int nsize,
                 const int nnzB1, //const double B1[],
                 const int B1rowId[], const int B1colSum[],
                 const int nnzB2, //const double B2[], 
                 const int B2rowId[], const int B2colSum[]);

int sparseUpdateAndCheckPattern(const double tol, const int nsize,
                                const int nnzA, const double A[],
                                const int ArowId[], const int AcolSum[],
                                int *nnzB, double B[], int BrowId[],
                                int BcolSum[], const double gamma,
                                const int threshType, bool doPatternCheck,
                                bool strictSamePattern);

int sparse_jac_v(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy,
                 void *user_data, N_Vector tmp);

int jac_full_dense(long int N, realtype t, N_Vector y, N_Vector fy,
			DlsMat Jac, void *user_data,
			N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int jac_full_sparse(int N, realtype t, N_Vector y, N_Vector fy,
      double user_jac[], void *user_data,
      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int jac_full_sparse_divided_diff(int N, realtype t, N_Vector y, N_Vector fy,
                    double user_jac[], void *user_data,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
