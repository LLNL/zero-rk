#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "lapack_wrapper.h"
#include "polynomial.h"

int GetPolynomialRoots(const int degree,
                       const double coefficients[],
                       double real_roots[],
                       double imag_roots[])
{
  const int matrix_size = degree*degree;
  int lapack_info=0;
  int j,k,m;
  double *A=NULL;

  if(degree < 1) {
    printf("ERROR: In GetPolynomialRoots(...),\n");
    printf("       no roots computed for a degree %d polynomial.\n",degree);
    fflush(stdout);
    return NO_LAPACK_CALL;
  }

  if(fabs(coefficients[degree])<MIN_COEFFICIENT) {
    printf("ERROR: In GetPolynomialRoots(...),\n");
    printf("       largest degree coefficient c[%d]*x^%d = %.18g\n",
           degree,degree,coefficients[degree]);
    printf("       will produce ill-conditioned system.\n");
    printf("       No roots computed.\n");
    fflush(stdout);
    return NO_LAPACK_CALL;

  }

  if(degree == 1) {
    // p(x) = c[0] + c[1]*x = 0  =>   x = -c[0]/c[1]
    real_roots[0] = -coefficients[0]/coefficients[1];
    imag_roots[0] = 0.0;
    return 0;
  }

  A = (double *)malloc(sizeof(double)*matrix_size);
  if(A == NULL) {
    printf("ERROR: In GetPolynomialRoots(...),\n");
    printf("       no roots computed.\n");
    printf("       Failed to allocate matrix size %d (bytes)\n",
           (int)sizeof(double)*matrix_size);
    printf("       (%d by %d)  sizeof(double) elements.\n",
           degree,degree);
    fflush(stdout);
    return NO_LAPACK_CALL;
  }
  // build the matrix
  m=0; // current matrix position
  for(j=0; j<(degree-1); ++j) { // column j
    for(k=0; k<degree; ++k) { // row k

      // lower sub-diagonal is equal to 1
      if(j+1 == k) {
        A[m] = 1.0;
      } else {
        A[m] = 0.0;
      }
      ++m; 
    } // for row-k
  }  // for column-j

  // Note that the last column is the negative of the monic polynomial
  // coefficients. The monic polynomial with the same roots as the target
  // polynomial is simply normalized by the largest degree coefficient
  for(k=0; k<degree; ++k) { // row k
    A[m] = -coefficients[k]/coefficients[degree];
    ++m;
  }
  //for(j=0; j<degree; ++j) { // column j
  //  for(k=0; k<degree; ++k) { // row k
  //    printf("A[%d][%d]: %10.4f\n",k,j,A[j*degree+k]);
  //  }
  //}
  
  lapack_info = LapackDGEEV('N',        // don't compute  left eigenvectors
                            'N',        // don't compute right eigenvectors
                            degree,     // # matrix rows
                            A,          // matrix to compute eigenvalues
                            degree,     // lda leading dimension of A
                            real_roots, // real eigenvalue part
                            imag_roots, // imaginary eigenvalue part
                            NULL,       // matrix VL for left eigenvectors
                            degree,     // ldvl leading dimension of VL
                            NULL,       // matrix VR for right eigenvectors
                            degree);    // ldvr leading dimension of VR

  if(A != NULL) {
    free(A);
  }
  return lapack_info;
  
}

int GetPolynomialExtrema(const int degree,
                         const double coefficients[],
                         const double imag_root_rtol,
                         const double imag_root_atol,
                         double x_extrema[],
                         double p_extrema[])
{
  int num_extrema=0;
  int error_flag;
  int true_degree;
  int j;
  double *real_roots;
  double *imag_roots;
  double *derivative;

  true_degree = -1;
  for(j=degree; j>=0; --j) {
    if(fabs(coefficients[j]) > MIN_COEFFICIENT) {
      true_degree  = j;
      break;
    }
  }
  if(true_degree == -1) {
    printf("WARNING: In GetPolynomialExtrema(...),\n");
    printf("         zero polynomial given.\n");
    fflush(stdout);
    return 0;
  }
  if(true_degree < 2) {
    // no local extrema possible for linear or constant polynomial
    return 0;  
  }
  
  derivative = (double *)malloc(sizeof(double)*true_degree);
  real_roots = (double *)malloc(sizeof(double)*(true_degree-1));
  imag_roots = (double *)malloc(sizeof(double)*(true_degree-1));

  // compute the derivative of the polynomial
  for(j=0; j<true_degree; ++j) {
    derivative[j] = coefficients[j+1]*(double)(j+1);
  }
  
  error_flag = GetPolynomialRoots(true_degree-1,
                                  derivative,
                                  real_roots,
                                  imag_roots);


  if(error_flag != 0) {
    printf("ERROR: In GetPolynomialExtrema(...),\n");
    printf("       GetPolynomialRoots(...) returned error flag = %d.\n",
           error_flag);
    printf("       Could not determine roots for the derivative of the\n");
    printf("       polynomial p(x).\n\n");
    printf("          p(x) = c[0] + c[1]*x + c[2]*x^2 + ...\n\n");
    for(j=0; j<=degree; ++j) {
      printf("            c[%d] = %.18g\n",j,coefficients[j]);
    }
    printf("       with derivative p'(x) = d[0] + d[1]*x + ...\n\n");
    for(j=0; j<true_degree; ++j) {
      printf("            d[%d] = %.18g\n",j,derivative[j]);
    }
    if(derivative != NULL) {
      free(derivative);
    }
    if(real_roots != NULL) {
      free(real_roots);
    }
    if(imag_roots != NULL) {
      free(imag_roots);
    }
    fflush(stdout);
    return -1;
  }
  // search for any real roots in the derivative indicating a local extrema
  for(j=0; j<(true_degree-1); ++j) {

    if(fabs(imag_roots[j]) <= fabs(real_roots[j]*imag_root_rtol) ||
       fabs(imag_roots[j]) <= imag_root_atol) {

      x_extrema[num_extrema] = real_roots[j];
      p_extrema[num_extrema] = EvaluatePolynomial(true_degree,
                                                  coefficients,
                                                  real_roots[j]);
      ++num_extrema;
    }
  }
  if(derivative != NULL) {
    free(derivative);
  }
  if(real_roots != NULL) {
    free(real_roots);
  }
  if(imag_roots != NULL) {
    free(imag_roots);
  }
  fflush(stdout);

  return num_extrema;
}

double EvaluatePolynomial(const int degree,
                          const double coefficients[],
                          const double x)
{
  int j;
  double polynomial_sum;

  if(degree < 0) {
    return 0.0;
  }

  polynomial_sum= coefficients[degree];
  for(j=degree-1; j>=0; --j) {
    polynomial_sum *= x;
    polynomial_sum += coefficients[j]; 
  }
  return polynomial_sum;
}
