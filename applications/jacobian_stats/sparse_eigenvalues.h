#ifndef SPARSE_EIGENVALUES_H_
#define SPARSE_EIGENVALUES_H_


int GetEigenvalues(const int num_rows,
                   const int num_nonzeros,
                   const double sparse_matrix[],
                   const int row_id[], 
                   const int column_sum[],
                   double real_eigenvalues[],
                   double imag_eigenvalues[]);


#endif

