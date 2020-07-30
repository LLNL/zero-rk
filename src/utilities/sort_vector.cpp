#include <stdlib.h>

#include <limits>

#include "sort_vector.h"
namespace zerork {
namespace utilities {

int CompareSortableVector(const void *aptr, const void *bptr) {
  SortableVector a = *(SortableVector *)aptr;
  SortableVector b = *(SortableVector *)bptr;
  double a_max = SortableVectorMax(&a);
  double b_max = SortableVectorMax(&b);
  if(a_max < b_max) {
    return -1;
  } else if(a_max > b_max) {
    return 1;
  }
  return 0;
}
int CompareSortableVectorDescending(const void *aptr, const void *bptr)
{return -CompareSortableVector(aptr,bptr);}


void SortByVectorMax(const int num_vectors,
                     SortableVector *list)
{
  qsort(&list[0],
        num_vectors,
        sizeof(list[0]),
        CompareSortableVectorDescending);
        

} 

// Return the largest vector element stored in the SortableVector struct. If
// the member SortableVector.v_double (type std::vector<double>) has size of
// zero, then the minimum finite double precision value is returned.
double SortableVectorMax(const SortableVector *sort_vec)
{
  const int length = sort_vec->v_double.size();
  if(length >= 1) {
    double sort_vec_max = sort_vec->v_double[0];
    for(int j=1; j<length; ++j) {
      if(sort_vec_max < sort_vec->v_double[j]) {
         sort_vec_max = sort_vec->v_double[j];
      }
    }
    return sort_vec_max;
  } else {
    return std::numeric_limits<double>::min();
  }
}

} // end namespace utilities
} // end namespace zerork
