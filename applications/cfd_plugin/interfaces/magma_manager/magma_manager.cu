#include "magma_manager.h"
#include "../../cuda_err_check.h"


template<typename T>
magma_manager<T>::magma_manager() :
  n_(-1),
  num_batches_(-1),
  factored_(false)
{
  magma_init();
  int device_id;
  magma_getdevice(&device_id);
  magma_queue_create(device_id, &magma_queue_);
}

template<typename T>
magma_manager<T>::~magma_manager()
{
  if(factored_) {
    FreeDeviceMemory();
  }
  magma_queue_destroy(magma_queue_);
  magma_finalize();
}

template<typename T>
void magma_manager<T>::setup_memory()
{
  if(factored_) {
    FreeDeviceMemory();
  }
  AllocateDeviceMemory();
}


template<typename T>
void magma_manager<T>::AllocateDeviceMemory()
{
  cudaDeviceSynchronize();
  cuda_err_check(cudaGetLastError());

  cuda_err_check(cudaMalloc((void**)&matrix_inverse_dev_,sizeof(T)*(n_*n_*num_batches_)));
  cuda_err_check(cudaMalloc((void**)&matrix_inverse_pointers_dev_,sizeof(T*)*num_batches_));
  cuda_err_check(cudaMalloc((void**)&matrix_pointers_dev_,sizeof(T*)*num_batches_));
  cuda_err_check(cudaMalloc((void**)&info_dev_,sizeof(int)*num_batches_));
  cuda_err_check(cudaMalloc((void**)&ipiv_dev_,sizeof(int)*n_*num_batches_));
  cuda_err_check(cudaMalloc((void**)&tmp_dev_,sizeof(T)*num_batches_*n_));
  cuda_err_check(cudaMalloc((void**)&tmp_pointers_dev_,sizeof(T*)*num_batches_));
  cuda_err_check(cudaMalloc((void**)&ipiv_pointers_dev_,sizeof(int*)*num_batches_));

  data_ptrs_.resize(num_batches_);

  std::vector<int*> tmpi_ptrs(num_batches_);
  std::vector<T*> tmp_ptrs(num_batches_);

  for(int j = 0; j < num_batches_; ++j) {
    data_ptrs_[j] = matrix_inverse_dev_ + j*n_*n_;
  }
  cudaMemcpy(matrix_inverse_pointers_dev_, data_ptrs_.data(), sizeof(T*)*num_batches_, cudaMemcpyHostToDevice);
  cuda_err_check(cudaGetLastError());

  for(int j = 0; j < num_batches_; ++j) {
    tmp_ptrs[j] = tmp_dev_ + j*n_;
  }
  cudaMemcpy(tmp_pointers_dev_, tmp_ptrs.data(), sizeof(T*)*num_batches_, cudaMemcpyHostToDevice);
  cuda_err_check(cudaGetLastError());

  for(int j = 0; j < num_batches_; ++j) {
    tmpi_ptrs[j] = ipiv_dev_ + j*n_;
  }
  cudaMemcpy(ipiv_pointers_dev_, tmpi_ptrs.data(), sizeof(int*)*num_batches_, cudaMemcpyHostToDevice);
  cuda_err_check(cudaGetLastError());
}

template<typename T>
void magma_manager<T>::FreeDeviceMemory()
{
  cudaFree(matrix_inverse_dev_);
  cudaFree(matrix_inverse_pointers_dev_);
  cudaFree(matrix_pointers_dev_);
  cudaFree(info_dev_);
  cudaFree(ipiv_dev_);
  cudaFree(ipiv_pointers_dev_);
  cudaFree(tmp_dev_);
  cudaFree(tmp_pointers_dev_);
}

template<>
void magma_manager<double>::getrf_batched() {
  magma_dgetrf_batched(n_, /* number of rows per block */
                       n_, /* number of columns per block */
                       matrix_pointers_dev_,
                       n_, /* leading dimension of each block */
                       ipiv_pointers_dev_,
                       info_dev_,
                       num_batches_,
                       magma_queue_);
}

template<>
void magma_manager<cuDoubleComplex>::getrf_batched() {
  magma_zgetrf_batched(n_, /* number of rows per block */
                       n_, /* number of columns per block */
                       matrix_pointers_dev_,
                       n_, /* leading dimension of each block */
                       ipiv_pointers_dev_,
                       info_dev_,
                       num_batches_,
                       magma_queue_);
}

template<>
void magma_manager<double>::getri_batched() {
  magma_dgetri_outofplace_batched(n_, /* order of block */
                                  matrix_pointers_dev_,
                                  n_, /* leading dimension of each block */
                                  ipiv_pointers_dev_,
                                  matrix_inverse_pointers_dev_,
                                  n_, /* leading dimension of each block of inverse */
                                  info_dev_,
                                  num_batches_,
                                  magma_queue_);
}

template<>
void magma_manager<cuDoubleComplex>::getri_batched() {
  magma_zgetri_outofplace_batched(n_, /* order of block */
                                  matrix_pointers_dev_,
                                  n_, /* leading dimension of each block */
                                  ipiv_pointers_dev_,
                                  matrix_inverse_pointers_dev_,
                                  n_, /* leading dimension of each block of inverse */
                                  info_dev_,
                                  num_batches_,
                                  magma_queue_);
}


template<typename T>
int magma_manager<T>::factor_invert(int num_batches, int n, T* values) {
  if(n != n_ || num_batches != num_batches_) {
    n_ = n;
    num_batches_ = num_batches;
    setup_memory();
  }
  if(values == NULL) {
    return 1;
  }

  bool need_tx = false;
  for(int j = 0; j < num_batches_; ++j) {
    if(data_ptrs_[j] != values + j*n_*n_) {
      data_ptrs_[j] = values + j*n_*n_;
      need_tx = true;
    }
  }
  if(need_tx) {
    cudaMemcpy(matrix_pointers_dev_, data_ptrs_.data(), sizeof(T*)*num_batches_, cudaMemcpyHostToDevice);
  }

  getrf_batched();
  getri_batched();


  int ierr = 0;
#ifdef ZERORK_FULL_DEBUG
  info_.resize(num_batches_);
  cuda_err_check(cudaMemcpy(info_.data(), info_dev_, num_batches_*sizeof(int), cudaMemcpyDeviceToHost));
  //Check for errors
  // factor_error > 0, singular matrix, zero diagonal at row,col = factor_error
  // factor_error = 0, success
  // factor_error < 0, illegal input
  for(int i=0; i < num_batches_; ++i) {
    if(info_[i]!=0) {
      ierr = info_[i];
      break;
    }
  }
#endif

  factored_ = true;
  return ierr;
}

template<typename T>
int magma_manager<T>::factor_lu(int num_batches, int n, T* values) {
  if(n != n_ || num_batches != num_batches_) {
    n_ = n;
    num_batches_ = num_batches;
    setup_memory();
  }
  if(values == NULL) {
    return 1;
  }

  bool need_tx = false;
  for(int j = 0; j < num_batches_; ++j) {
    if(data_ptrs_[j] != values + j*n_*n_) {
      data_ptrs_[j] = values + j*n_*n_;
      need_tx = true;
    }
  }
  if(need_tx) {
    cudaMemcpy(matrix_pointers_dev_, data_ptrs_.data(), sizeof(T*)*num_batches_, cudaMemcpyHostToDevice);
  }

  getrf_batched();

  int ierr = 0;
#ifdef ZERORK_FULL_DEBUG
  info_.resize(num_batches_);
  cuda_err_check(cudaMemcpy(info_.data(), info_dev_, num_batches_*sizeof(int), cudaMemcpyDeviceToHost));
  //Check for errors
  // factor_error > 0, singular matrix, zero diagonal at row,col = factor_error
  // factor_error = 0, success
  // factor_error < 0, illegal input
  for(int i=0; i < num_batches_; ++i) {
    if(info_[i]!=0) {
      ierr = info_[i];
      break;
    }
  }
#endif

  factored_ = true;
  return ierr;
}


//The following modified from cuda sdk-5.0
#define TRANSPOSE_TILE_DIM    32
#define TRANSPOSE_BLOCK_ROWS  8

template<typename T>
static __global__ void MAGMA_MANAGER_TransposeNoBankConflicts(T *odata, const T *idata, const int width, const int height)
{
    __shared__ T tile[TRANSPOSE_TILE_DIM][TRANSPOSE_TILE_DIM+1];
    int xIndex,yIndex,index_in,index_out;

    xIndex = blockIdx.x * TRANSPOSE_TILE_DIM + threadIdx.x;
    yIndex = blockIdx.y * TRANSPOSE_TILE_DIM + threadIdx.y;
    index_in = xIndex + (yIndex)*width;

    for (int i=0; i<TRANSPOSE_TILE_DIM; i+=TRANSPOSE_BLOCK_ROWS)
    {
        if(xIndex < width && yIndex+i < height){
        tile[threadIdx.y+i][threadIdx.x] = idata[index_in+i*width];}
    }

    __syncthreads();

    xIndex = blockIdx.y * TRANSPOSE_TILE_DIM + threadIdx.x;
    yIndex = blockIdx.x * TRANSPOSE_TILE_DIM + threadIdx.y;
    index_out = xIndex + (yIndex)*height;

    for (int i=0; i<TRANSPOSE_TILE_DIM; i+=TRANSPOSE_BLOCK_ROWS)
    {
        if(yIndex+i < width && xIndex < height){
        odata[index_out+i*height] = tile[threadIdx.x][threadIdx.y+i];}
    }
}

template<typename T>
void magma_manager<T>::cuda_transpose(T* odata, const T* idata, const int width, const int height)
{
    // Put df/dy in "normal" order
    dim3 nBlocks2D,nThreads2D;
    nThreads2D.x = TRANSPOSE_TILE_DIM;
    nThreads2D.y = TRANSPOSE_BLOCK_ROWS;
    nBlocks2D.x = (width+TRANSPOSE_TILE_DIM-1)/TRANSPOSE_TILE_DIM;
    nBlocks2D.y = (height+TRANSPOSE_TILE_DIM-1)/TRANSPOSE_TILE_DIM;
    MAGMA_MANAGER_TransposeNoBankConflicts<T><<<nBlocks2D,nThreads2D>>>(odata,idata,width,height);
#ifdef ZERORK_FULL_DEBUG
    cuda_err_check( cudaPeekAtLastError() );
    cuda_err_check( cudaDeviceSynchronize() );
#endif
}


namespace {
template<typename T>
void __global__ MAGMA_MANAGER_cuda_bdmv_kernel
(
    const int mtx_block_size,
    const int num_mtx_blocks,
    const T* A_dev,
    const T* X_dev ,
    T * Y_dev
)
{
  int tidx = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = gridDim.x*blockDim.x;
  for( ; tidx < num_mtx_blocks*mtx_block_size; tidx += stride)
  {
    int local_row   = tidx % mtx_block_size;
    int local_block = tidx / mtx_block_size;
    T Y_dev_accum = 0.0;
    for(int i = 0; i < mtx_block_size; ++i) //columns
    {
      int data_idx = mtx_block_size*mtx_block_size*local_block + mtx_block_size*i + local_row;
      Y_dev_accum += A_dev[data_idx]*X_dev[i+local_block*mtx_block_size];
    }
    Y_dev[local_row+local_block*mtx_block_size] = Y_dev_accum;
  }
}

template<>
void __global__ MAGMA_MANAGER_cuda_bdmv_kernel
(
    const int mtx_block_size,
    const int num_mtx_blocks,
    const cuDoubleComplex* A_dev,
    const cuDoubleComplex* X_dev ,
    cuDoubleComplex * Y_dev
)
{
  int tidx = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = gridDim.x*blockDim.x;
  for( ; tidx < num_mtx_blocks*mtx_block_size; tidx += stride)
  {
    int local_row   = tidx % mtx_block_size;
    int local_block = tidx / mtx_block_size;
    cuDoubleComplex Y_dev_accum = make_cuDoubleComplex(0.0,0.0);
    for(int i = 0; i < mtx_block_size; ++i) //columns
    {
      int data_idx = mtx_block_size*mtx_block_size*local_block + mtx_block_size*i + local_row;
      //Y_dev_accum += A_dev[data_idx]*X_dev[i+local_block*mtx_block_size];
      Y_dev_accum = cuCadd(Y_dev_accum, cuCmul(A_dev[data_idx],X_dev[i+local_block*mtx_block_size]));
    }
    Y_dev[local_row+local_block*mtx_block_size] = Y_dev_accum;
  }
}
} //anonymous namespace

template<typename T>
int magma_manager<T>::cuda_bdmv(int n, int nbatch, T* A_dev, T* B_dev, T* Y_dev)
{
  int threads = std::min(n*nbatch,1024);
  int blocks=(nbatch*n+threads-1)/threads;
  MAGMA_MANAGER_cuda_bdmv_kernel<T><<<blocks,threads>>>(n, nbatch, A_dev, B_dev, Y_dev);
#ifdef ZERORK_FULL_DEBUG
  cuda_err_check(cudaPeekAtLastError());
  cuda_err_check(cudaDeviceSynchronize());
#endif
  return 0;  
}

template<typename T>
int magma_manager<T>::solve_invert(int num_batches, int n, const T* rhs, T* soln) {
  if(n != n_ || num_batches != num_batches_) {
    return 1;
  }

  // Transpose rhs into soln
  cuda_transpose(soln,rhs,num_batches_,n_);

  // Block-diagonal matrix vector multiplication
  cuda_bdmv(n_, num_batches_, matrix_inverse_dev_, soln, tmp_dev_);

  // Put tmp back into block order
  cuda_transpose(soln,tmp_dev_,n_,num_batches_);

  return(0);
}

template<>
void magma_manager<double>::getrs_batched() {
  // Magma forward and back substitution
  magma_dgetrs_batched(MagmaNoTrans,
                       n_, /* order of the matrix */
                       1, /* number of right hand sides */
                       matrix_pointers_dev_,
                       n_, /* leading dimension of A */
                       ipiv_pointers_dev_,
                       tmp_pointers_dev_, /* right hand side (input), solution (output) */
                       n_, /* leading dimension of b */
                       num_batches_,
                       magma_queue_);
}

template<>
void magma_manager<cuDoubleComplex>::getrs_batched() {
  // Magma forward and back substitution
  magma_zgetrs_batched(MagmaNoTrans,
                       n_, /* order of the matrix */
                       1, /* number of right hand sides */
                       matrix_pointers_dev_,
                       n_, /* leading dimension of A */
                       ipiv_pointers_dev_,
                       tmp_pointers_dev_, /* right hand side (input), solution (output) */
                       n_, /* leading dimension of b */
                       num_batches_,
                       magma_queue_);
}

template<typename T>
int magma_manager<T>::solve_lu(int num_batches, int n, const T* rhs, T* soln) {
  if(n != n_ || num_batches != num_batches_) {
    return 1;
  }

  // Transpose rhs into tmp_dev_
  cuda_transpose(tmp_dev_,rhs,num_batches_,n_);

  getrs_batched();

  // Put tmp back into block order
  cuda_transpose(soln,tmp_dev_,n_,num_batches_);

  return(0);
}

template class magma_manager<double>;
template class magma_manager<cuDoubleComplex>;

