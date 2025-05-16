
#include <cuComplex.h>
#include <stdio.h>

#include "cuda_la_manager.h"

template<typename T>
cuda_la_manager<T>::cuda_la_manager() 
 :
   use_lu_(false)
{
};

template<typename T>
int cuda_la_manager<T>::factor(int num_batches, int n, T* values) {
 if(use_lu_) {
   return(this->factor_lu(num_batches, n, values));
 } else {
   return(this->factor_invert(num_batches, n, values));
 }
}

template<typename T>
int cuda_la_manager<T>::solve(int num_batches, int n, const T* rhs, T* soln) {
 if(use_lu_) {
   return(this->solve_lu(num_batches, n, rhs, soln));
 } else {
   return(this->solve_invert(num_batches, n, rhs, soln));
 }
}


template class cuda_la_manager<double>;
template class cuda_la_manager<cuDoubleComplex>;

