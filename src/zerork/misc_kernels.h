#ifndef ZERORK_MISC_KERNELS_H
#define ZERORK_MISC_KERNELS_H

namespace zerork {


void setDoubleArrayVal
(
    double *array,
    const double val,
    const int idxStart = 0,
    const int idxEnd = 1,
    const cudaStream_t stream = 0
);

void printCudaArray(const int n, const double *A_dev);


}// namespace zerork

#endif
