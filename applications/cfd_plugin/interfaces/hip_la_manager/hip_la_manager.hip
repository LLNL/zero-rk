
#include <stdio.h>
#include <hip/hip_complex.h>
#include "hip_la_manager.h"


template<typename T>
hip_la_manager<T>::hip_la_manager() 
 :
   use_lu_(false)
{
};

template<typename T>
int hip_la_manager<T>::factor(int num_batches, int n, T* values) {
 if(use_lu_) {
   return(this->factor_lu(num_batches, n, values));
 } else {
   return(this->factor_invert(num_batches, n, values));
 }
}

template<typename T>
int hip_la_manager<T>::solve(int num_batches, int n, const T* rhs, T* soln) {
 if(use_lu_) {
   return(this->solve_lu(num_batches, n, rhs, soln));
 } else {
   return(this->solve_invert(num_batches, n, rhs, soln));
 }
}


template class hip_la_manager<double>;
template class hip_la_manager<hipDoubleComplex>;

