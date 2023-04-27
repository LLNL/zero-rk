#include "hip/hip_runtime.h"
#include "hipblas_manager.h"

#include "../../gpu_err_check.h"


hipblas_manager::hipblas_manager() :
  n_(-1),
  num_batches_(-1),
  factored_(false)
{
  hipblasCreate(&hipblas_handle_);
}

hipblas_manager::~hipblas_manager()
{
  if(factored_) {
    FreeDeviceMemory();
    hipblasDestroy(hipblas_handle_);
  }
}

void hipblas_manager::setup_memory()
{
  if(factored_) {
    FreeDeviceMemory();
  }
  AllocateDeviceMemory();
}


void hipblas_manager::AllocateDeviceMemory()
{
  hipDeviceSynchronize();
  gpu_err_check(hipGetLastError());

  gpu_err_check(hipMalloc((void**)&matrix_inverse_dev_,sizeof(double)*(n_*n_*num_batches_)));
  gpu_err_check(hipMalloc((void**)&matrix_inverse_pointers_dev_,sizeof(double*)*num_batches_));
  gpu_err_check(hipMalloc((void**)&matrix_pointers_dev_,sizeof(double*)*num_batches_));
  gpu_err_check(hipMalloc((void**)&info_dev_,sizeof(int)*num_batches_));
  gpu_err_check(hipMalloc((void**)&tmp_dev_,sizeof(double)*num_batches_*n_));
  gpu_err_check(hipMalloc((void**)&tmp_pointers_dev_,sizeof(double*)*num_batches_));

  data_ptrs_.resize(num_batches_);
  tmp_ptrs_.resize(num_batches_);

  for(int j = 0; j < num_batches_; ++j) {
    data_ptrs_[j] = matrix_inverse_dev_ + j*n_*n_;
  }
  hipMemcpy(matrix_inverse_pointers_dev_, data_ptrs_.data(), sizeof(double*)*num_batches_, hipMemcpyHostToDevice);
  gpu_err_check(hipGetLastError());

  for(int j = 0; j < num_batches_; ++j) {
    tmp_ptrs_[j] = tmp_dev_ + j*n_;
  }
  hipMemcpy(tmp_pointers_dev_, tmp_ptrs_.data(), sizeof(double*)*num_batches_, hipMemcpyHostToDevice);
  gpu_err_check(hipGetLastError());
}

void hipblas_manager::FreeDeviceMemory()
{
  hipFree(matrix_inverse_dev_);
  hipFree(matrix_inverse_pointers_dev_);
  hipFree(matrix_pointers_dev_);
  hipFree(info_dev_);
  hipFree(tmp_dev_);
  hipFree(tmp_pointers_dev_);
}

int hipblas_manager::factor_invert(int num_batches, int n, double* values) {
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
    hipMemcpy(matrix_pointers_dev_, data_ptrs_.data(), sizeof(double*)*num_batches_, hipMemcpyHostToDevice);
  }

  int lda = n_;
  int* ipiv = NULL; //Turns off pivoting
  hipblasDgetrfBatched(hipblas_handle_, n_,
                      matrix_pointers_dev_, lda,
                      ipiv, info_dev_, num_batches_);

  int ldc = n_;
  double* const* const_matrix_pointers_dev = (double* const*) matrix_pointers_dev_;
  hipblasDgetriBatched(hipblas_handle_, n_, const_matrix_pointers_dev,
                      lda, ipiv, matrix_inverse_pointers_dev_,
                      ldc, info_dev_, num_batches_);

  int ierr = 0;
#ifdef ZERORK_FULL_DEBUG
  info_.resize(num_batches_);
  gpu_err_check(hipMemcpy(info_.data(), info_dev_, num_batches_*sizeof(int), hipMemcpyDeviceToHost));
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

int hipblas_manager::factor_lu(int num_batches, int n, double* values) {
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
    hipMemcpy(matrix_pointers_dev_, data_ptrs_.data(), sizeof(double*)*num_batches_, hipMemcpyHostToDevice);
  }

  int lda = n_;
  int* ipiv = NULL; //Turns off pivoting
  hipblasDgetrfBatched(hipblas_handle_, n_,
                      matrix_pointers_dev_, lda,
                      ipiv, info_dev_, num_batches_);

  int ierr = 0;
#ifdef ZERORK_FULL_DEBUG
  info_.resize(num_batches_);
  gpu_err_check(hipMemcpy(info_.data(), info_dev_, num_batches_*sizeof(int), hipMemcpyDeviceToHost));
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

static __global__ void CUBLAS_MANAGER_TransposeNoBankConflicts(double *odata, const double *idata, const int width, const int height)
{
    __shared__ double tile[TRANSPOSE_TILE_DIM][TRANSPOSE_TILE_DIM+1];
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

void hipblas_manager::gpu_transpose(double* odata, const double* idata, const int width, const int height)
{
    // Put df/dy in "normal" order
    dim3 nBlocks2D,nThreads2D;
    nThreads2D.x = TRANSPOSE_TILE_DIM;
    nThreads2D.y = TRANSPOSE_BLOCK_ROWS;
    nBlocks2D.x = (width+TRANSPOSE_TILE_DIM-1)/TRANSPOSE_TILE_DIM;
    nBlocks2D.y = (height+TRANSPOSE_TILE_DIM-1)/TRANSPOSE_TILE_DIM;
    CUBLAS_MANAGER_TransposeNoBankConflicts<<<nBlocks2D,nThreads2D>>>(odata,idata,width,height);
#ifdef ZERORK_FULL_DEBUG
    gpu_err_check( hipPeekAtLastError() );
    gpu_err_check( hipDeviceSynchronize() );
#endif
}


static void __global__ CUBLAS_MANAGER_gpu_bdmv_kernel
(
    const int mtx_block_size,
    const int num_mtx_blocks,
    const double* A_dev,
    const double* X_dev ,
    double * Y_dev
)
{
  int tidx = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = gridDim.x*blockDim.x;
  for( ; tidx < num_mtx_blocks*mtx_block_size; tidx += stride)
  {
    int local_row   = tidx % mtx_block_size;
    int local_block = tidx / mtx_block_size;
    double Y_dev_accum = 0.0;
    for(int i = 0; i < mtx_block_size; ++i) //columns
    {
      int data_idx = mtx_block_size*mtx_block_size*local_block + mtx_block_size*i + local_row;
      Y_dev_accum += A_dev[data_idx]*X_dev[i+local_block*mtx_block_size];
    }
    Y_dev[local_row+local_block*mtx_block_size] = Y_dev_accum;
  }
}

int hipblas_manager::gpu_bdmv(int n, int nbatch, double* A_dev, double* B_dev, double* Y_dev)
{
  int threads = std::min(n*nbatch,1024);
  int blocks=(nbatch*n+threads-1)/threads;
  CUBLAS_MANAGER_gpu_bdmv_kernel<<<blocks,threads>>>(n, nbatch, A_dev, B_dev, Y_dev);
#ifdef ZERORK_FULL_DEBUG
  gpu_err_check(hipPeekAtLastError());
  gpu_err_check(hipDeviceSynchronize());
#endif
  return 0;  
}

int hipblas_manager::solve_invert(int num_batches, int n, const double* rhs, double* soln) {
  if(n != n_ || num_batches != num_batches_) {
    return 1;
  }

  // Transpose rhs into soln
  gpu_transpose(soln,rhs,num_batches_,n_);

  // Block-diagonal matrix vector multiplication
  gpu_bdmv(n_, num_batches_, matrix_inverse_dev_, soln, tmp_dev_);

  // Put tmp back into block order
  gpu_transpose(soln,tmp_dev_,n_,num_batches_);

  return(0);
}

int hipblas_manager::solve_lu(int num_batches, int n, const double* rhs, double* soln) {
  if(n != n_ || num_batches != num_batches_) {
    return 1;
  }

  // Transpose rhs into tmp_dev_
  gpu_transpose(tmp_dev_,rhs,num_batches_,n_);

  // CUBLAS forward and back substitution
  int* ipiv = NULL; //Turns off pivoting
  int lda = n_;
  int ldb = n_;
  int info = 0;
  hipblasDgetrsBatched(hipblas_handle_, HIPBLAS_OP_N, n_, 1,
                      matrix_pointers_dev_, lda,
                      ipiv, tmp_pointers_dev_, ldb, &info, num_batches_);

  // Put tmp back into block order
  gpu_transpose(soln,tmp_dev_,n_,num_batches_);

  return(0);
}

