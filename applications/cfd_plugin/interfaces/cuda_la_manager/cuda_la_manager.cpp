
#include <stdio.h>
#include "cuda_la_manager.h"

cuda_la_manager::cuda_la_manager() 
 :
   use_lu_(false)
{
};

int cuda_la_manager::factor(int num_batches, int n, double* values) {
 if(use_lu_) {
   return(this->factor_lu(num_batches, n, values));
 } else {
   return(this->factor_invert(num_batches, n, values));
 }
}

int cuda_la_manager::solve(int num_batches, int n, const double* rhs, double* soln) {
 if(use_lu_) {
   return(this->solve_lu(num_batches, n, rhs, soln));
 } else {
   return(this->solve_invert(num_batches, n, rhs, soln));
 }
}


