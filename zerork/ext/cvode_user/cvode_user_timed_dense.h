/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006/11/29 00:05:06 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Header file for the CVODE dense linear solver CVLAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _CVUSERTIMEDLAPACK_H
#define _CVUSERTIMEDLAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

int CVUserTimedLapackDense(void *cvode_mem, int N);
int CVUserTimedGetSetupTime(void *cvode_mem, double * setup_time);
int CVUserTimedGetSolveTime(void *cvode_mem, double * solve_time);
int CVUserTimedGetJacEvalTime(void *cvode_mem, double * jac_eval_time);

#ifdef __cplusplus
}
#endif

#endif
