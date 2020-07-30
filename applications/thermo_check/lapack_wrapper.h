#ifndef LAPACK_WRAPPER_H_
#define LAPACK_WRAPPER_H_

#ifdef __cplusplus
extern "C" {
#endif

// DGELS solves overdetermined or underdetermined real linear systems
// involving an M-by-N matrix A, or its transpose, using a QR or LQ
// factorization of A.  It is assumed that A has full rank.
//
// The following options are provided:
//
// 1. If TRANS = 'N' and m >= n:  find the least squares solution of
//    an overdetermined system, i.e., solve the least squares problem
//                 minimize || B - A*X ||.
//
// 2. If TRANS = 'N' and m < n:  find the minimum norm solution of
//    an underdetermined system A * X = B.
//
// 3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
//    an undetermined system A**T * X = B.
//
// 4. If TRANS = 'T' and m < n:  find the least squares solution of
//    an overdetermined system, i.e., solve the least squares problem
//                 minimize || B - A**T * X ||.
//
// Several right hand side vectors b and solution vectors x can be
// handled in a single call; they are stored as the columns of the
// M-by-NRHS right hand side matrix B and the N-by-NRHS solution
// matrix X.

int LapackDGELS(char transpose,
                int num_rows,
                int num_cols,
                int num_rhs,
                double A[],
                int lda,
                double B[],
                int ldb);

// DGGLSE solves the linear equality-constrained least squares (LSE)
// problem:
//
//         minimize || c - A*x ||_2   subject to   B*x = d
//
// where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
// M-vector, and d is a given P-vector. It is assumed that
// P <= N <= M+P, and
//
//          rank(B) = P and  rank( (A) ) = N.
//                               ( (B) )
//
// These conditions ensure that the LSE problem has a unique solution,
// which is obtained using a generalized RQ factorization of the
// matrices (B, A) given by
//
//     B = (0 R)*Q,   A = Z*T*Q.

int LapackDGGLSE(int num_rows_A,
                 int num_cols_A,
                 int num_rows_B,
                 double A[],
                 int lda,
                 double B[],
                 int ldb,
                 double c[],
                 double d[],
                 double solution_x[]);
// DGEEV computes for an N-by-N real nonsymmetric matrix A, the
// eigenvalues and, optionally, the left and/or right eigenvectors.
//
// The right eigenvector v(j) of A satisfies
//                  A * v(j) = lambda(j) * v(j)
// where lambda(j) is its eigenvalue.
// The left eigenvector u(j) of A satisfies
//               u(j)**H * A = lambda(j) * u(j)**H
// where u(j)**H denotes the conjugate-transpose of u(j).
//
// The computed eigenvectors are normalized to have Euclidean norm
// equal to 1 and largest component real.

int LapackDGEEV(char find_left_eigenvectors,
                char find_right_eigenvectors,
                int num_rows,
                double A[],
                int lda,
                double real_eigenvalues[],
                double imag_eigenvalues[],
                double left_eigenvectors[],
                int ldvl,
                double right_eigenvectors[],
                int ldvr);


#ifdef __cplusplus
}
#endif

#endif
