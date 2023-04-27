
#include <stdio.h>
#include "hip_la_manager.h"

hip_la_manager::hip_la_manager() 
 :
   use_lu_(false)
{
};

int hip_la_manager::factor(int num_batches, int n, double* values) {
 if(use_lu_) {
   return(this->factor_lu(num_batches, n, values));
 } else {
   return(this->factor_invert(num_batches, n, values));
 }
}

int hip_la_manager::solve(int num_batches, int n, const double* rhs, double* soln) {
 if(use_lu_) {
   return(this->solve_lu(num_batches, n, rhs, soln));
 } else {
   return(this->solve_invert(num_batches, n, rhs, soln));
 }
}


