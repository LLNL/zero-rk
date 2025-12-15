/* -----------------------------------------------------------------
 * Zero-RK Modification of Sundials Magma Linear Solver to include 
 * re-ordering of equations
 *
 * -----------------------------------------------------------------
 * Original Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 * ----------------------------------------------------------------*/

#ifndef _ZRKLINSOL_MAGMADENSE_H
#define _ZRKLINSOL_MAGMADENSE_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_memory.h>
#include <sundials/sundials_nvector.h>

#if defined(SUNDIALS_MAGMA_BACKENDS_CUDA)
#define HAVE_CUBLAS
#elif defined(SUNDIALS_MAGMA_BACKENDS_HIP)
#define HAVE_HIP
#endif
#include <magma_v2.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------
 * MAGMA dense implementation of SUNLinearSolver
 * ----------------------------------------------- */

struct _ZRKLinearSolverContent_MagmaDense {
  int             last_flag;
  booleantype     async;
  sunindextype    N;
  SUNMemory       pivots;
  SUNMemory       pivotsarr;
  SUNMemory       dpivotsarr;
  SUNMemory       infoarr;
  SUNMemory       rhsarr;
  SUNMemoryHelper memhelp;
  magma_queue_t   q;
  SUNMatrix       Ainv;
  booleantype     use_lu;
};

typedef struct _ZRKLinearSolverContent_MagmaDense *ZRKLinearSolverContent_MagmaDense;


SUNDIALS_EXPORT SUNLinearSolver ZRKLinSol_MagmaDense(N_Vector y, SUNMatrix A);

SUNDIALS_EXPORT int ZRKLinSol_MagmaDense_SetAsync(SUNLinearSolver S, booleantype onoff);

SUNDIALS_EXPORT SUNLinearSolver_Type ZRKLinSolGetType_MagmaDense(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID ZRKLinSolGetID_MagmaDense(SUNLinearSolver S);
SUNDIALS_EXPORT int ZRKLinSolInitialize_MagmaDense(SUNLinearSolver S);
SUNDIALS_EXPORT int ZRKLinSolSetup_MagmaDense(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int ZRKLinSolSolve_MagmaDense(SUNLinearSolver S, SUNMatrix A,
                                              N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT sunindextype ZRKLinSolLastFlag_MagmaDense(SUNLinearSolver S);
SUNDIALS_EXPORT int ZRKLinSolSpace_MagmaDense(SUNLinearSolver S,
                                              long int *lenrwLS,
                                              long int *leniwLS);
SUNDIALS_EXPORT int ZRKLinSolFree_MagmaDense(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
