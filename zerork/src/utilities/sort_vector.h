#ifndef SORT_VECTOR_H_
#define SORT_VECTOR_H_

#include <vector>

namespace zerork {

typedef struct {
  int id;
  std::vector<double> v_double;
} SortableVector;

int CompareSortableVector(const void *aptr, const void *bptr);

int CompareSortableVectorDescending(const void *aptr, const void *bptr);

void SortByVectorMax(const int num_vectors,
                     SortableVector *list); 

double SortableVectorMax(const SortableVector *A);

} // end of namespace zerork
#endif
