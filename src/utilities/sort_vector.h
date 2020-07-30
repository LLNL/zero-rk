#ifndef SORT_VECTOR_H_
#define SORT_VECTOR_H_

#include <vector>

namespace zerork {
namespace utilities {

typedef struct {
  int id;
  std::vector<double> v_double;
} SortableVector;

int CompareSortableVector(const void *aptr, const void *bptr);

int CompareSortableVectorDescending(const void *aptr, const void *bptr);

void SortByVectorMax(const int num_vectors,
                     SortableVector *list); 

double SortableVectorMax(const SortableVector *A);

} // end namespace utilities
} // end of namespace zerork
#endif
