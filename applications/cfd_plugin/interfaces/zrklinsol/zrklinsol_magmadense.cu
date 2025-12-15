/* -----------------------------------------------------------------
 * Zero-RK Modification of Sundials Magma Linear Solver to include
 * re-ordering of equations
 *
 * Original Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 * ----------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <algorithm> // for std::min

#include "zrklinsol_magmadense.h"
#include <sunmatrix/sunmatrix_magmadense.h>
#include <sundials/sundials_math.h>

/* Interfaces to match 'realtype' with the correct MAGMA functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define xgetrf magma_dgetrf_gpu
#define xgetrf_batched magma_dgetrf_batched
#define xgetrs magma_dgetrs_gpu
#define xgetrs_batched magma_dgetrs_batched
#define xgetri_batched magma_dgetri_outofplace_batched
#define xset_pointer magma_dset_pointer
#define xgeam magmablas_dgeam
#define xtranspose magmablas_dtranspose
#define xcopy magma_dcopy
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define xgetrf magma_sgetrf_gpu
#define xgetrf_batched magma_sgetrf_batched
#define xgetrs magma_sgetrs_gpu
#define xgetrs_batched magma_sgetrs_batched
#define xgetri_batched magma_sgetri_outofplace_batched
#define xset_pointer magma_sset_pointer
#define xgeam magmablas_sgeam
#define xtranspose magmablas_stranspose
#define xcopy magma_scopy
#else
#error  Incompatible realtype for MAGMA
#endif

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * MAGMADENSE solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define MAGMADENSE_CONTENT(S) ( (ZRKLinearSolverContent_MagmaDense)(S->content) )
#define MHELP(S)              ( MAGMADENSE_CONTENT(S)->memhelp )
#define QUEUE(S)              ( MAGMADENSE_CONTENT(S)->q )
#define PIVOTS(S)             ( (sunindextype*)MAGMADENSE_CONTENT(S)->pivots->ptr )
#define PIVOTSARRAY(S)        ( (sunindextype**)MAGMADENSE_CONTENT(S)->pivotsarr->ptr )
#define RHSARRAY(S)           ( (realtype**)MAGMADENSE_CONTENT(S)->rhsarr->ptr )
#define INFOARRAY(S)          ( (sunindextype*)MAGMADENSE_CONTENT(S)->infoarr->ptr )
#define LASTFLAG(S)           ( MAGMADENSE_CONTENT(S)->last_flag )
#define ASYNCHRONOUS(S)       ( MAGMADENSE_CONTENT(S)->async)

/*
 * ----------------------------------------------------------------------------
 * Implementation specific routines
 * ----------------------------------------------------------------------------
 */

/*
 * Constructor functions
 */

SUNLinearSolver ZRKLinSol_MagmaDense(N_Vector y, SUNMatrix Amat)
{
  int retval = 0;
  SUNLinearSolver S;
  ZRKLinearSolverContent_MagmaDense content;
  SUNMatrixContent_MagmaDense A;
  sunindextype M, nblocks;

  /* Check inputs */
  if (y == NULL || Amat == NULL)
    return(NULL);

  if (y->ops == NULL || Amat->ops == NULL)
    return(NULL);

  if (y->ops->nvgetlength == NULL || y->ops->nvgetdevicearraypointer == NULL ||
      Amat->ops->getid == NULL)
    return(NULL);

  /* Check compatibility with supplied SUNMatrix */
  if (SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE)
    return(NULL);

  if (Amat->content == NULL)
    return(NULL);

  A = (SUNMatrixContent_MagmaDense) Amat->content;

  /* Check that the matrix is square */
  if (A->M != A->N)
    return(NULL);

  M = A->M;
  nblocks = A->nblocks;

  /* Check that the matirx and vector dimensions agree */
  if (M*nblocks != N_VGetLength(y))
    return(NULL);

  /* Create the linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty();
  if (S == NULL) return(NULL);

  /* Attach operations */
  S->ops->gettype    = ZRKLinSolGetType_MagmaDense;
  S->ops->getid      = ZRKLinSolGetID_MagmaDense;
  S->ops->initialize = ZRKLinSolInitialize_MagmaDense;
  S->ops->setup      = ZRKLinSolSetup_MagmaDense;
  S->ops->solve      = ZRKLinSolSolve_MagmaDense;
  S->ops->lastflag   = ZRKLinSolLastFlag_MagmaDense;
  S->ops->space      = ZRKLinSolSpace_MagmaDense;
  S->ops->free       = ZRKLinSolFree_MagmaDense;

  /* Create content */
  content = NULL;
  content = (ZRKLinearSolverContent_MagmaDense) malloc(sizeof(*content));
  if (content == NULL) { SUNLinSolFree(S); return(NULL); }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->last_flag = 0;
  content->async     = SUNTRUE;
  content->N         = M;
  content->pivots    = NULL;
  content->pivotsarr = NULL;
  content->infoarr   = NULL;
  content->rhsarr    = NULL;
  content->memhelp   = A->memhelp;
  content->q         = A->q;
  content->Ainv      = SUNMatClone(Amat);

  //Set boolean for LU/inverse solve
  content->use_lu = SUNFALSE;  // default
  if(getenv("ZERORK_REACTOR_USE_LU") != NULL) {
    int use_lu_env = atoi(getenv("ZERORK_REACTOR_USE_LU"));
    if(use_lu_env != 0) {
        content->use_lu = SUNTRUE;
    }
  }

  /* Allocate data */

  /* The pivots need to be in host memory when calling the non-batched methods,
     but in device memory for the batched methods. */
  retval = SUNMemoryHelper_Alloc(content->memhelp, &content->pivots,
                                 M*nblocks*sizeof(sunindextype),
                                 nblocks > 1 ? SUNMEMTYPE_DEVICE : SUNMEMTYPE_HOST);
  if (retval) { SUNLinSolFree(S); return(NULL); }

  /* If we have multiple blocks, then we need to allocate some extra
     pointer arrays needed when calling MAGMA batched methods. */
  if (nblocks > 1)
  {
    retval = SUNMemoryHelper_Alloc(content->memhelp, &content->pivotsarr,
                                   nblocks*sizeof(sunindextype*), SUNMEMTYPE_DEVICE);
    if (retval) { SUNLinSolFree(S); return(NULL); }

    /* Set the pivots array on the device */
    magma_iset_pointer((sunindextype**)content->pivotsarr->ptr, /* 2D output array */
                       (sunindextype*)content->pivots->ptr,  /* 1D input array */
                       1, /* leading dimension of input */
                       0, /* row */
                       0, /* column */
                       M, /* rows in a block */
                       nblocks, /* number of blocks */
                       content->q);

    /* We use pinned memory for the info array because we are going to
       check its values on the host and we need it to have fast transfers. */
    retval = SUNMemoryHelper_Alloc(content->memhelp, &content->infoarr,
                                   nblocks*sizeof(sunindextype), SUNMEMTYPE_PINNED);
    if (retval) { SUNLinSolFree(S); return(NULL); }

    retval = SUNMemoryHelper_Alloc(content->memhelp, &content->rhsarr,
                                   nblocks*sizeof(realtype*), SUNMEMTYPE_DEVICE);
    if (retval) { SUNLinSolFree(S); return(NULL); }
  }

  return(S);
}

/*
 * Set functions
 */

int ZRKLinSol_MagmaDense_SetAsync(SUNLinearSolver S, booleantype onoff)
{
  if (S == NULL) return SUNLS_MEM_NULL;
  ASYNCHRONOUS(S) = onoff;
  return SUNLS_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * Implementation of generic SUNLinearSolver operations.
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type ZRKLinSolGetType_MagmaDense(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}

SUNLinearSolver_ID ZRKLinSolGetID_MagmaDense(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_CUSTOM);
}

int ZRKLinSolInitialize_MagmaDense(SUNLinearSolver S)
{
  /* All solver-specific memory has already been allocated */
  if (S == NULL) return SUNLS_MEM_NULL;
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(SUNLS_SUCCESS);
}

int ZRKLinSolSetup_MagmaDense(SUNLinearSolver S, SUNMatrix A)
{
  /* Check for valid inputs */
  if (S == NULL) return SUNLS_MEM_NULL;

  if (A == NULL)
  {
    LASTFLAG(S) = SUNLS_MEM_NULL;
    return(SUNLS_MEM_NULL);
  }

  /* Ensure that A is a magma dense matrix */
  if (SUNMatGetID(A) != SUNMATRIX_MAGMADENSE)
  {
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return(SUNLS_ILL_INPUT);
  }

  sunindextype ier = 0;
  sunindextype M = SUNMatrix_MagmaDense_BlockRows(A);
  sunindextype nblocks = SUNMatrix_MagmaDense_NumBlocks(A);

  SUNMatrix Ainv = MAGMADENSE_CONTENT(S)->Ainv;

  /* Call MAGMA to do LU factorization of A */
  if (nblocks > 1)
  {
#ifndef SUNDIALS_MAGMA_USE_GETRF_LOOP
    xgetrf_batched(M, /* number of rows per block */
                   M, /* number of columns per block */
                   SUNMatrix_MagmaDense_BlockData(A),
                   M, /* leading dimension of each block */
                   PIVOTSARRAY(S),
                   INFOARRAY(S),
                   nblocks,
                   QUEUE(S));
    if(!MAGMADENSE_CONTENT(S)->use_lu){
      //Only solve for the direct inverse
      xgetri_batched(M, /* order of block */
                     SUNMatrix_MagmaDense_BlockData(A),
                     M, /* leading dimension of each block */
                     PIVOTSARRAY(S),
                     SUNMatrix_MagmaDense_BlockData(Ainv),
                     M, /* leading dimension of each block of inverse */
                     INFOARRAY(S),
                     nblocks,
                     QUEUE(S));
    }
#else
    realtype** blocks = SUNMatrix_MagmaDense_BlockData(A);
    for (int k = 0; k < nblocks; k++)
    {
      xgetrf(M, /* number of rows */
             M, /* number of columns */
             blocks[k],
             M, /* leading dimension of A */
             PIVOTSARRAY(S)[k],
             &INFOARRAY(S)[k]);
    }
#endif

    if (!ASYNCHRONOUS(S))
    {
      magma_queue_sync(QUEUE(S));
      /* Check if there were any failures when factoring */
      for (sunindextype k = 0; k < nblocks; k++)
      {
        if (INFOARRAY(S)[k] < 0) ier = INFOARRAY(S)[k];
        if (INFOARRAY(S)[k] > 0)
        {
          ier = INFOARRAY(S)[k];
          break;
        }
      }
    }
  }
  else
  {
    xgetrf(M, /* number of rows */
           M, /* number of columns */
           SUNMatrix_MagmaDense_Data(A),
           M, /* leading dimension of A */
           PIVOTS(S),
           &ier);
    if (!ASYNCHRONOUS(S)) magma_queue_sync(QUEUE(S));
  }

  LASTFLAG(S) = ier;
  if (ier > 0) return(SUNLS_LUFACT_FAIL);
  if (ier < 0) return(SUNLS_PACKAGE_FAIL_UNREC);
  return(SUNLS_SUCCESS);
}

namespace {
static void __global__ ZRKLINSOL_magma_bdmv_kernel
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
}

int ZRKLinSolSolve_MagmaDense(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                              N_Vector b, realtype tol)
{
  /* Check for valid inputs */
  if (S == NULL) return(SUNLS_MEM_NULL);

  if ( (A == NULL) || (x == NULL) || (b == NULL) )
  {
    LASTFLAG(S) = SUNLS_MEM_NULL;
    return(SUNLS_MEM_NULL);
  }

  /* Ensure that A is a magma dense matrix */
  if (SUNMatGetID(A) != SUNMATRIX_MAGMADENSE)
  {
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return(SUNLS_ILL_INPUT);
  }

  int ier = 0;
  sunindextype M = SUNMatrix_MagmaDense_BlockRows(A);
  sunindextype nblocks = SUNMatrix_MagmaDense_NumBlocks(A);

  SUNMatrix Ainv = MAGMADENSE_CONTENT(S)->Ainv;

  /* Access x data array */
  realtype* xdata = N_VGetDeviceArrayPointer(x);
  if (xdata == NULL)
  {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(SUNLS_MEM_FAIL);
  }

  /* Access x data array */
  realtype* bdata = N_VGetDeviceArrayPointer(b);
  if (bdata == NULL)
  {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(SUNLS_MEM_FAIL);
  }

  /* Call MAGMA to solve the linear system */
  if (nblocks > 1)
  {
    if(MAGMADENSE_CONTENT(S)->use_lu) {
      //Need to transpose tog et the order right for xset_pointer.
      xtranspose(nblocks, M, bdata, nblocks, xdata, M, QUEUE(S));

      /* Set up the pointers for each of the blocks */
      xset_pointer(RHSARRAY(S),
                   xdata, /* 1D input array */
                   1, /* leading dimension of input */
                   0, /* rows */
                   0, /* cols */
                   M, /* number of rows in block */
                   nblocks,
                   QUEUE(S));

      /* Now, solve the batch system */
      xgetrs_batched(MagmaNoTrans,
                     M, /* order of the matrix */
                     1, /* number of right hand sides */
                     SUNMatrix_MagmaDense_BlockData(A),
                     M, /* leading dimension of A */
                     PIVOTSARRAY(S),
                     RHSARRAY(S), /* right hand side (input), solution (output) */
                     M, /* leading dimension of b */
                     nblocks,
                     QUEUE(S));

      if (!ASYNCHRONOUS(S)) magma_queue_sync(QUEUE(S));

      //Need to transpose back to what ZeroRK expects
      xtranspose(M, nblocks, xdata, M, bdata, nblocks, QUEUE(S));

      N_VScale(ONE, b, x);

    } else {
      //Use the direct solve
      xtranspose(nblocks, M, bdata, nblocks, xdata, M, QUEUE(S));
      int threads = std::min(M * nblocks, 1024);
      int blocks = (nblocks * M + threads - 1) / threads;
      ZRKLINSOL_magma_bdmv_kernel<<<blocks,threads>>>(M, nblocks,
              SUNMatrix_MagmaDense_Data(Ainv), xdata, bdata
      );
      xtranspose(M, nblocks, bdata, M, xdata, nblocks, QUEUE(S));
      if (!ASYNCHRONOUS(S)) magma_queue_sync(QUEUE(S));
    }
  }
  else
  {
    xgetrs(MagmaNoTrans,
           M, /* order of the matrix */
           1, /* number of right hand sides */
           SUNMatrix_MagmaDense_Data(A),
           M, /* leading dimension of A */
           PIVOTS(S),
           xdata, /* right hand side (input), solution (output) */
           M, /* leading dimension of x */
           &ier);
  }
  if(!ASYNCHRONOUS(S)) magma_queue_sync(QUEUE(S));

  LASTFLAG(S) = ier;
  return((ier < 0) ? SUNLS_PACKAGE_FAIL_UNREC : SUNLS_SUCCESS);
}

sunindextype ZRKLinSolLastFlag_MagmaDense(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  if (S == NULL) return(-1);
  return(LASTFLAG(S));
}

int ZRKLinSolSpace_MagmaDense(SUNLinearSolver S,
                              long int *lenrwLS,
                              long int *leniwLS)
{
  *lenrwLS = 0;
  *leniwLS = 2 + MAGMADENSE_CONTENT(S)->N;
  return(SUNLS_SUCCESS);
}

int ZRKLinSolFree_MagmaDense(SUNLinearSolver S)
{
  /* return if S is already free */
  if (S == NULL) return(SUNLS_SUCCESS);

  /* delete items from contents, then delete generic structure */
  if (S->content)
  {
    if (MAGMADENSE_CONTENT(S)->pivots)
      SUNMemoryHelper_Dealloc(MHELP(S), MAGMADENSE_CONTENT(S)->pivots);
    if (MAGMADENSE_CONTENT(S)->pivotsarr)
      SUNMemoryHelper_Dealloc(MHELP(S), MAGMADENSE_CONTENT(S)->pivotsarr);
    if (MAGMADENSE_CONTENT(S)->infoarr)
      SUNMemoryHelper_Dealloc(MHELP(S), MAGMADENSE_CONTENT(S)->infoarr);
    if (MAGMADENSE_CONTENT(S)->rhsarr)
      SUNMemoryHelper_Dealloc(MHELP(S), MAGMADENSE_CONTENT(S)->rhsarr);
    SUNMatDestroy(MAGMADENSE_CONTENT(S)->Ainv);
  }
  if (S->ops)
  {
    free(S->ops);
    S->ops = NULL;
  }
  free(S);
  S = NULL;
  return(SUNLS_SUCCESS);
}
