/*
 * -----------------------------------------------------------------
 * Zero-RK Modification of Sundials Dense Matrix Module to use
 * BLAS routines
 *
 * Original Programmer(s): Daniel Reynolds @ SMU
 *                         David Gardner @ LLNL
 * Based on code sundials_direct.h by: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 * -----------------------------------------------------------------
 */

#ifndef _ZRKLAPACKMATRIX_DENSE_H
#define _ZRKLAPACKMATRIX_DENSE_H

#include <stdio.h>
#include <sundials/sundials_matrix.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* ----------------------------------
 * Dense implementation of SUNMatrix
 * ---------------------------------- */
  
struct _SUNMatrixContent_ZRKLapackDense {
  sunindextype M;
  sunindextype N;
  realtype *data;
  sunindextype ldata;
  realtype **cols;
};

typedef struct _SUNMatrixContent_ZRKLapackDense *SUNMatrixContent_ZRKLapackDense;

/* ------------------------------------
 * Macros for access to SUNMATRIX_DENSE
 * ------------------------------------ */

#define ZM_CONTENT_D(A)     ( (SUNMatrixContent_ZRKLapackDense)(A->content) )

#define ZM_ROWS_D(A)        ( ZM_CONTENT_D(A)->M )

#define ZM_COLUMNS_D(A)     ( ZM_CONTENT_D(A)->N )

#define ZM_LDATA_D(A)       ( ZM_CONTENT_D(A)->ldata )

#define ZM_DATA_D(A)        ( ZM_CONTENT_D(A)->data )

#define ZM_COLS_D(A)        ( ZM_CONTENT_D(A)->cols )

#define ZM_COLUMN_D(A,j)    ( (ZM_CONTENT_D(A)->cols)[j] )

#define ZM_ELEMENT_D(A,i,j) ( (ZM_CONTENT_D(A)->cols)[j][i] )

/* ---------------------------------------
 * Exported Functions for SUNMATRIX_DENSE
 * --------------------------------------- */

SUNDIALS_EXPORT SUNMatrix ZRKDenseLapackMatrix(sunindextype M, sunindextype N);

SUNDIALS_EXPORT void ZRKDenseLapackMatrix_Print(SUNMatrix A, FILE* outfile);

SUNDIALS_EXPORT sunindextype ZRKDenseLapackMatrix_Rows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype ZRKDenseLapackMatrix_Columns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype ZRKDenseLapackMatrix_LData(SUNMatrix A);
SUNDIALS_EXPORT realtype* ZRKDenseLapackMatrix_Data(SUNMatrix A);
SUNDIALS_EXPORT realtype** ZRKDenseLapackMatrix_Cols(SUNMatrix A);
SUNDIALS_EXPORT realtype* ZRKDenseLapackMatrix_Column(SUNMatrix A, sunindextype j);

SUNDIALS_EXPORT SUNMatrix_ID SUNMatGetID_ZRKLapackDense(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix SUNMatClone_ZRKLapackDense(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_ZRKLapackDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatZero_ZRKLapackDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatCopy_ZRKLapackDense(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAdd_ZRKLapackDense(realtype c, SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAddI_ZRKLapackDense(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvec_ZRKLapackDense(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int SUNMatSpace_ZRKLapackDense(SUNMatrix A, long int *lenrw, long int *leniw);

  
#ifdef __cplusplus
}
#endif

#endif
