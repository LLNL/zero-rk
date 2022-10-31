#ifndef ZERORK_CUDA_DEFS_H
#define ZERORK_CUDA_DEFS_H

#include <cuda_runtime.h>

#define MAX_NREACS 5
#define MAX_NPRODS 5

#define WARP_SIZE 32
#define WARPS_PER_BLOCK 16
#define THREADS_PER_BLOCK (WARP_SIZE*WARPS_PER_BLOCK)
#define MAX_THREADS_PER_BLOCK 1024

namespace zerork {

void checkCudaError(cudaError_t err, const char *msg = NULL);

} // namespace zerork
#endif
