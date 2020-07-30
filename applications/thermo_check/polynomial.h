#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#ifdef __cplusplus
extern "C" {
#endif

const int NO_LAPACK_CALL = -100000000;
const double MIN_COEFFICIENT = 1.0e-300;

// Compute the roots of a polynomial p(x) of degree d, with real-valued 
// coefficients c[0], c[1], ..., c[d] such that
//
// p(x) = c[0] + c[1]*x + c[2]*x^2 + ... + c[d]*x^d (eq 1)
//
// The roots are computed by solving the eigenvalue problem for the companion
// matrix that has a characteristic polynomial with the same roots as p(x).
//
// Inputs:
//   degree           polynomial degree (one or greater)
//   coefficients[i]  polynomial coefficient c[i] in (eq 1),
//                    array length must be degree+1
// Outputs:
//   real_roots[i]    real part of the polynomial roots,
//                    array length must be degree
//   imag_roots[i]    imaginary part of the polynomial roots,
//                    array length must be degree
//
// Return: NO_LAPACK_CALL if the lapack eigenvalue function DGEEV could not be
//         called, otherwise the lapack info code returned from the 
//         eigenvalue calculation DGEEV. The lapack info codes are as follows:
//
//         info = 0:  successful exit
//                    < 0:  if INFO = -i, the i-th argument had an 
//                          illegal value.
//                    > 0:  if INFO = i, the QR algorithm failed to compute
//                          all the eigenvalues, and no eigenvectors have 
//                          been computed; elements i+1:N of WR and WI contain 
//                          eigenvalues which have converged.
	    
int GetPolynomialRoots(const int degree,
                       const double coefficients[],
                       double real_roots[],
                       double imag_roots[]);

int GetPolynomialExtrema(const int degree,
                         const double coefficients[],
                         const double imag_root_rtol,
                         const double imag_root_atol,
                         double x_extrema[],
                         double p_extrema[]);
double EvaluatePolynomial(const int degree,
                          const double coefficients[],
                          const double x);

#ifdef __cplusplus
}
#endif

#endif
