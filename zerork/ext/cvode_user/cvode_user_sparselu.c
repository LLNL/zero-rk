#include "cvode_user_sparselu.h"

// CVUserSuperLU sets various function fields and initializes the memory
// associated with the linear solver. Also needs to set setupNonNull in
// *cvode_mem to TRUE.
int CVUserSuperLU(void  *cvode_mem,
                  const int num_eq,
                  const int num_nz,
                  const int row_id[],
                  const int col_start[],
                  int (*jac_func_ptr)(int,realtype, N_Vector, N_Vector,
                                      double *, void *, N_Vector, N_Vector,
                                      N_Vector))
{
  int j;
  int flag;
  int *sparse_diagonal_id;
  CVodeMem cv_mem;
  CVUserSuperLUMem cv_user_lmem;

  //printf("In CVUserSuperLU\n"); fflush(stdout);

  if(cvode_mem == NULL) {
    // TODO: issue cvode memory error
    printf("CVode Memory pointer is NULL\n"); fflush(stdout);
    return -1;
  }
  cv_mem = (CVodeMem) cvode_mem;

  // Free any previously attached linear solver memory
  if(cv_mem->cv_lfree != NULL) {
    cv_mem->cv_lfree(cv_mem);
  }
 
  // set the four main linear solver functions
  cv_mem->cv_linit  = &LInitUserSuperLU;
  cv_mem->cv_lsetup = &LSetupUserSuperLU;
  cv_mem->cv_lsolve = &LSolveUserSuperLU;
  cv_mem->cv_lfree  = &LFreeUserSuperLU;

  cv_mem->cv_setupNonNull = TRUE;

  // set the memory for CVUserSuperLUMemRec
  cv_user_lmem = NULL;
  cv_user_lmem = (CVUserSuperLUMem) malloc(sizeof(struct CVUserSuperLUMemRec));
  if(cv_user_lmem == NULL) {
    // TODO: issue linear solver memory allocation error
    printf("User linear solver memory pointer is NULL\n"); fflush(stdout);
    return -1;
  }
  // check that the sparse Jacobian has a complete diagonal
  sparse_diagonal_id = (int *)malloc(sizeof(int)*num_eq);
  flag = HasCompleteDiagonal(num_eq,
                             num_nz,
                             row_id,
                             col_start,
                             sparse_diagonal_id);

  if(flag != 1) {
    // TODO: issue sparse jacobian diagonal 
    printf("Sparse Jacobian does not have a diagonal element in row %d\n",
           1-flag);
    fflush(stdout);
    return flag-1;
  } 
  // set the user-supplied Jacobian function
  cv_user_lmem->UserJacobianFunction = jac_func_ptr;

  // initialize Jacobian related data
  // :
  // :
  cv_user_lmem->num_equations = num_eq;
  cv_user_lmem->num_nonzeros  = num_nz;
  cv_user_lmem->num_jac       = 0; // number of jacobians calculated
  cv_user_lmem->num_factor    = 0; // number of factorizations
  cv_user_lmem->num_solve     = 0; // number of backward substitutions
  

  // allocate memory for the sparse Jacobian matrix user data
  cv_user_lmem->row_index       = (int *)malloc(sizeof(int)*num_nz);
  cv_user_lmem->column_sum      = (int *)malloc(sizeof(int)*(num_eq+1));
  cv_user_lmem->sparse_identity = (double *)malloc(sizeof(double)*num_nz);
  cv_user_lmem->J_user          = (double *)malloc(sizeof(double)*num_nz);
  cv_user_lmem->M_user          = (double *)malloc(sizeof(double)*num_nz);
  cv_user_lmem->b_user          = (double *)malloc(sizeof(double)*num_eq);
  cv_user_lmem->x_user          = (double *)malloc(sizeof(double)*num_eq);

  // copy the sparsity pattern information
  memcpy(cv_user_lmem->row_index,
         row_id,
         num_nz*sizeof(int));
  memcpy(cv_user_lmem->column_sum,
         col_start,
         (num_eq+1)*sizeof(int));

  // create a sparese identity matrix
  for(j=0; j<num_nz; ++j) {
    cv_user_lmem->sparse_identity[j]=0.0;
  }
  for(j=0; j<num_eq; ++j) {
    cv_user_lmem->sparse_identity[sparse_diagonal_id[j]] = 1.0;
  } 

  // SuperLU specific initializations
  cv_user_lmem->row_permutation = (int *)malloc(sizeof(int)*num_eq);
  cv_user_lmem->col_permutation = (int *)malloc(sizeof(int)*num_eq);
  cv_user_lmem->col_elim_tree   = (int *)malloc(sizeof(int)*num_eq);
  cv_user_lmem->Rvec_internal   = (double *)malloc(sizeof(double)*num_eq);
  cv_user_lmem->Cvec_internal   = (double *)malloc(sizeof(double)*num_eq);

  /* Set the default input options:
	options.Fact = DOFACT;
        options.Equil = YES;
    	options.ColPerm = COLAMD;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
	options.DiagPivotThresh = 1.0;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES;
  */
  set_default_options(&(cv_user_lmem->optionSLU));
  // MMD_AT_PLUS_A performs best out of stock SuperLU options for iso-octane
  // preconditioner.
  //cv_user_lmem->optionSLU.ColPerm = MMD_AT_PLUS_A;

  // create a compressed column matrix
  // SLU_NC ==
  // SLU_D  == data type double precision
  // SLU_GE == matrix structure general
  dCreate_CompCol_Matrix(&(cv_user_lmem->Mslu),
			 num_eq,              //cv_user_lmem->num_equations,
			 num_eq,              //cv_user_lmem->num_equations,
			 num_nz,              //cv_user_lmem->num_nonzeros,
			 cv_user_lmem->M_user,
			 cv_user_lmem->row_index,
			 cv_user_lmem->column_sum,
			 SLU_NC,SLU_D,SLU_GE);

  dCreate_Dense_Matrix(&(cv_user_lmem->Bslu),
                       num_eq,              //cv_user_lmem->num_equations,
                       1,
                       cv_user_lmem->b_user,
                       num_eq,              //cv_user_lmem->num_equations,
                       SLU_DN,SLU_D,SLU_GE);
  dCreate_Dense_Matrix(&(cv_user_lmem->Xslu),
                       num_eq,              //cv_user_lmem->num_equations,
                       1,
                       cv_user_lmem->x_user,
                       num_eq,              //cv_user_lmem->num_equations,
                       SLU_DN,SLU_D,SLU_GE);

  StatInit(&(cv_user_lmem->statSLU)); // initialize SuperLU statistics
  

  // attach the linear solver memory to the full integrator memory (CVodeMem)
  cv_mem->cv_lmem = cv_user_lmem;
  
  // free local memory
  free(sparse_diagonal_id);

  return 0;
}

// Four basic functions to provide an alternate linear solver
//  1. Initialization
int LInitUserSuperLU(CVodeMem cv_mem)
{
  CVUserSuperLUMem cv_user_lmem;

  //printf("In LInitUserSuperLU(...)\n"); fflush(stdout);

  cv_user_lmem = (CVUserSuperLUMem) cv_mem->cv_lmem;

  // clear counters
  cv_user_lmem->num_jac = 0;
  cv_user_lmem->num_factor = 0;
  cv_user_lmem->num_solve = 0;
  // In cvDenseInit the other tasks involved establishing the jacobian
  // function and associated data depending on whether or not the user provided
  // a jacobian function.  In it's current form, the user is expected to
  // provide the jacobian function.  If this is not the case in the future,
  // some function/data initialization/reset should be considered here.
  cv_user_lmem->last_flag = CV_SUCCESS;
  return 0;
}

//  2. Setup - prepares the linear solver for subsequent calls to LSolve 
// Following the steps in cvDenseSetup(...):
//
// This routine does the setup operations for the dense linear solver.
// It makes a decision whether or not to call the Jacobian evaluation
// routine based on various state variables, and if not it uses the 
// saved copy.  In any case, it constructs the Newton matrix 
// M = I - gamma*J, updates counters, and calls the dense LU 
// factorization routine.
int LSetupUserSuperLU(CVodeMem cv_mem,
                      int convfail,
                      N_Vector ypred,
                      N_Vector fpred,
                      booleantype *jcurPtr,
                      N_Vector vtemp1,
                      N_Vector vtemp2,
                      N_Vector vtemp3)
{
  // CVode controls defined in $(CVODE_ROOT)/src/cvode/cvode_direct_impl.h
  // These constants are not available through the installed headers.
  //   CVD_DGMAX - maximum change in gamma between Jacobian evaluations
  //   CVD_MSBJ  - maximum number of steps between Jacobian evaluations
  const double CVD_DGMAX=0.2; // [0.2 in cvode_direct_impl.h]
  const int CVD_MSBJ=50;      // [50  in cvode_direct_impl.h]
  
  booleantype jbad, jok;
  realtype dgamma;
  int factor_error;
  int retval;
  CVUserSuperLUMem cv_user_lmem = (CVUserSuperLUMem) cv_mem->cv_lmem;
  // for memcpy() use the following:
  //const size_t copy_size=cv_user_lmem->num_nonzeros*sizeof(double);
  // otherwise:
  const int copy_size = cv_user_lmem->num_nonzeros;
  //int lapack_size = cv_user_lmem->num_nonzeros;
  //int lapack_increment = 1;
  double lapack_scalar;
  // SuperLU dummy arguments and return arguments
  int lwork=0;
  void *work=NULL;
  double rpg,rcond;
  double ferr[1],berr[1]; // length is the number of RHS to be solved 
                          // simultaneously
  int j;                  

  //printf("In LSetupUserSuperLU(...)\n"); fflush(stdout);

  // Use nst, gamma/gammap and convfail to set the Jacobian evaluation flag jok
  // cv_mem->cv_nst is the number of internal integration steps.
  // cv_mem->gamma is the current value of the timestep divided by l1 (ell-one)
  //   which is used to scale the ODE Jacobian J (df/dy) for the nonlinear
  //   solver Jacobian M = I - gamma * J.
  // cv_mem->gammap is the previous value of cv_mem->gamma.
  dgamma = fabs(cv_mem->cv_gamma/cv_mem->cv_gammap) - 1.0;
  jbad = (cv_mem->cv_nst == 0) ||
         (cv_mem->cv_nst > cv_user_lmem->num_steps_last_jac  + CVD_MSBJ) ||
         ((convfail == CV_FAIL_BAD_J) && (dgamma < CVD_DGMAX)) ||
         (convfail == CV_FAIL_OTHER);
  jok = !jbad;
  //jok = FALSE;

  if(jok) {
    // If jok == TRUE, use a saved copy of J
    //memcpy(cv_user_lmem->M_user,
    //       cv_user_lmem->J_user,
    //       copy_size);
    //printf("jok == TRUE\n"); fflush(stdout);
    for(j=0; j<copy_size; ++j) {
      cv_user_lmem->M_user[j] = cv_user_lmem->J_user[j];
    }
    *jcurPtr = FALSE;
  } else {
    // If jok == FALSE, call the user specified sparse jacobian routine
    // for a new value of J
    //printf("jok == FALSE\n"); fflush(stdout);
    ++cv_user_lmem->num_jac; // number of jacobian evaluations
    cv_user_lmem->num_steps_last_jac = cv_mem->cv_nst;
    *jcurPtr = TRUE;

    // set the Jacobian to zero initially
    for(j=0; j<copy_size; ++j) {
      cv_user_lmem->M_user[j] = 0.0;
    }
 
    retval = (cv_user_lmem->UserJacobianFunction)(cv_user_lmem->num_equations,
                                                  cv_mem->cv_tn,
                                                  ypred,
                                                  fpred,
                                                  cv_user_lmem->M_user,
                                                  cv_mem->cv_user_data,
                                                  vtemp1,
                                                  vtemp2,
                                                  vtemp3);
    ++cv_user_lmem->num_jac;
    if(retval < 0) {
      // unrecoverable error returned from the user's jacobian function
      printf("ERROR: In LSetupUserSuperLU(...),\n");
      printf("       Unrecoverable error = %d returned from the user's\n",
             retval);
      printf("       jacobian function.\n");
      cv_user_lmem->last_flag = -5; // CVDLS_JACFUNC_UNRECVR;
      return -1;                    // same return as cvDenseSetup 
    }
    if(retval > 0) {
      // recoverable error returned from the user's jacobian function
      cv_user_lmem->last_flag = -6; // CVDLS_JACFUNC_RECVR;
      return 1;                     // same return as cvDenseSetup 
    }

    // if retval == 0 (CV_SUCCESS), then copy M_user to J_user for saving
    for(j=0; j<copy_size; ++j) {
      cv_user_lmem->J_user[j] = cv_user_lmem->M_user[j];
    }
      
  } // end of if(jok)-else block
  
  // At this point in the setup procedure, the M_user matrix contains the 
  // ODE system Jacobian J (which is stored in cv_user_lmem->J_user).
  // TODO: update these functions with the BLAS and LAPACK equivalents dscal
  //       and daxpy
  lapack_scalar = -cv_mem->cv_gamma;
  for(j=0; j<copy_size; ++j) {
    cv_user_lmem->M_user[j]*=lapack_scalar;
  }
  for(j=0; j<copy_size; ++j) {
    cv_user_lmem->M_user[j]+=cv_user_lmem->sparse_identity[j];
  }

  // perform LU factorization of M
  cv_user_lmem->optionSLU.Fact=SamePattern; // default except for first factor
  if(cv_user_lmem->num_factor > 0) {
    // You must destroy the L and U before each new factorization, or
    // you will have a memory leak
    Destroy_SuperNode_Matrix(&(cv_user_lmem->Lslu));
    Destroy_CompCol_Matrix(&(cv_user_lmem->Uslu));
  }
  else {
     // DOFACT indicates that the column & row permutation is called in
     // addition to the factorization
     // SamePattern indicates that the column permutation is reused
     // SamePattern_SameRowPerm indicates that the column and row permutations
     // are reused.  It is unclear if this option can be used reliably.
     cv_user_lmem->optionSLU.Fact=DOFACT;
     for(j=0; j<cv_user_lmem->num_equations; ++j) {
       cv_user_lmem->b_user[j] = 0.0;
     }
  }
  
  dgssvx(&(cv_user_lmem->optionSLU),
         &(cv_user_lmem->Mslu),
	 cv_user_lmem->col_permutation,
         cv_user_lmem->row_permutation,
	 cv_user_lmem->col_elim_tree,
         cv_user_lmem->equed,
         cv_user_lmem->Rvec_internal,
	 cv_user_lmem->Cvec_internal,
         &(cv_user_lmem->Lslu),
         &(cv_user_lmem->Uslu),
	 work,
         lwork,
         &(cv_user_lmem->Bslu),
         &(cv_user_lmem->Xslu),
         &rpg,
         &rcond,
	 ferr,
         berr,
         &(cv_user_lmem->mem_usage),
         &(cv_user_lmem->statSLU),
	 &factor_error);

  ++cv_user_lmem->num_factor;

  // factor_error > 0, singular matrix, zero diagonal at row,col = factor_error
  // factor_error = 0, success
  // factor_error < 0, illegal input
  cv_user_lmem->last_flag = factor_error;

  if(factor_error < 0) {
    printf("ERROR: In LSetupUserSuperLU(...),\n");
    printf("       SuperLU's dgssvx(...) returned error flag = %d\n",
           factor_error);
    printf("       Exiting Now.\n");
    exit(-1);
  }
  // Return 0 if LU was complete, otherwise return 1
  return factor_error;
}

// 3. Solve - solve the linear system
int LSolveUserSuperLU(CVodeMem cv_mem,
                      N_Vector b,
                      N_Vector weight,
                      N_Vector ycur,
                      N_Vector fcur)
{
  double *bptr = NV_DATA_S(b);
  CVUserSuperLUMem cv_user_lmem = (CVUserSuperLUMem) cv_mem->cv_lmem;
  // for memcpy() use the following:
  //const size_t copy_size=cv_user_lmem->num_nonzeros*sizeof(double);
  // otherwise:
  const int copy_size = cv_user_lmem->num_equations;
  // SuperLU dummy arguments and return arguments
  int lwork=0;
  void *work=NULL;
  double rpg,rcond;
  double ferr[1],berr[1]; // length is the number of RHS to be solved 
                          // simultaneously
  int solve_error;

  int j;

  //printf("In LSolveUserSuperLU(...)\n"); fflush(stdout);
  for(j=0; j<copy_size; ++j) {
    cv_user_lmem->b_user[j] = bptr[j];
  }

  cv_user_lmem->optionSLU.Fact=FACTORED;
  cv_user_lmem->Bslu.ncol=1; // in dlinsolx1.c example this is reset
                             // to the number of right hand sides
                             // not sure if this is still needed

  //backsolve with expert driver function
  dgssvx(&(cv_user_lmem->optionSLU),
         &(cv_user_lmem->Mslu),
	 cv_user_lmem->col_permutation,
         cv_user_lmem->row_permutation,
	 cv_user_lmem->col_elim_tree,
         cv_user_lmem->equed,
         cv_user_lmem->Rvec_internal,
	 cv_user_lmem->Cvec_internal,
         &(cv_user_lmem->Lslu),
         &(cv_user_lmem->Uslu),
	 work,
         lwork,
         &(cv_user_lmem->Bslu),
         &(cv_user_lmem->Xslu),
         &rpg,
         &rcond,
	 ferr,
         berr,
         &(cv_user_lmem->mem_usage),
         &(cv_user_lmem->statSLU),
	 &solve_error);

  ++cv_user_lmem->num_solve;

  // copy the solution from Xslu and x_user back to the b vector from the
  // function argument
  for(j=0; j<copy_size; ++j) {
    bptr[j] = cv_user_lmem->x_user[j];
  }
  // In cvDenseSolve the correction returned in b is scaled to account
  // for changes in gamma.  Chapter 7 in the CVode User's Guide does not
  // mention the need to scale the correction factor returned by lsolve
  // (see Section 7.3).
  // The following if-block is the same as found in cvDenseSolve:
  if ((cv_mem->cv_lmm == CV_BDF) && (cv_mem->cv_gamrat != 1.0)) {
    printf("INFO: In LSolveUserSuperLU(...),\n");
    printf("      applying gamma ratio correction to solution in b.\n");
    printf("      b = correction*b, correction = %.18g\n",
           (2.0/(1.0 + cv_mem->cv_gamrat)));
    fflush(stdout);
    N_VScale(2.0/(1.0 + cv_mem->cv_gamrat), b, b);
  }
  // the last flag and return values are the same as found in cvDenseSolve
  // factor_error > 0, singular matrix, zero diagonal at row,col = factor_error
  // factor_error = 0, success
  // factor_error < 0, illegal input
  cv_user_lmem->last_flag = solve_error;

  if(solve_error < 0) {
    printf("ERROR: In LSolveUserSuperLU(...),\n");
    printf("       SuperLU's dgssvx(...) returned error flag = %d\n",
           solve_error);
    printf("       Exiting Now.\n");
    exit(-1);
  }
  return solve_error;
}

// 4. Free any allocated memory
void LFreeUserSuperLU(CVodeMem cv_mem)
{
  CVUserSuperLUMem cv_user_lmem;

  //printf("In LFreeUserSuperLU(...)\n"); fflush(stdout);

  cv_user_lmem = (CVUserSuperLUMem) cv_mem->cv_lmem;

  free(cv_user_lmem->row_index);
  free(cv_user_lmem->column_sum);
  free(cv_user_lmem->sparse_identity);
  free(cv_user_lmem->J_user);
  free(cv_user_lmem->M_user);
  free(cv_user_lmem->b_user);
  free(cv_user_lmem->x_user);


  // SuperLU specific de-allocations
  free(cv_user_lmem->row_permutation);
  free(cv_user_lmem->col_permutation);
  free(cv_user_lmem->col_elim_tree);
  free(cv_user_lmem->Rvec_internal);
  free(cv_user_lmem->Cvec_internal);
  StatFree(&(cv_user_lmem->statSLU));
  Destroy_SuperMatrix_Store(&(cv_user_lmem->Bslu));
  Destroy_SuperMatrix_Store(&(cv_user_lmem->Xslu));

  if(cv_user_lmem->num_jac > 0) {
    Destroy_SuperNode_Matrix(&(cv_user_lmem->Lslu));
    Destroy_CompCol_Matrix(&(cv_user_lmem->Uslu));
  }
  Destroy_SuperMatrix_Store(&(cv_user_lmem->Mslu));

  // free struct CVUserSuperLUMemRec
  free(cv_user_lmem);
}


// Check the compressed column storage for all the diagonal elements.
// Return x: x == 1: all diagonal elements are represented
//           x <= 0: the (|x|+1)^th diagonal element was not found in the
//                   sparse matrix
// User can provide an array to store the diagonal terms location within
// the sparse data list via diag_sparse_id[].  If this is passed as a NULL
// pointer, then no data is stored.
int HasCompleteDiagonal(const int num_eq,
                        const int num_nz,
                        const int row_id[],
		        const int col_start[],
                        int diag_sparse_id[])
{
  int found_diag,sparse_id,num_in_col;
  int j,k;

  for(j=0; j<num_eq; ++j) {
    num_in_col = col_start[j+1]-col_start[j];
    found_diag = 0;
    for(k=0; k<num_in_col; ++k) {
      sparse_id = col_start[j]+k;
      if(row_id[sparse_id] == j) {
        found_diag = 1;
        // record the diagonal elements sparse index
        if(diag_sparse_id != NULL) {
          diag_sparse_id[j] = sparse_id;
        }
        break;
      }
    }
    if(found_diag == 0) {
      return -j;
    }  
  }
  return 1; // found all the diagonal elements
}

int CVUserSuperLUGetNumJacEvals(void *cv_mem,
                                int *num_jac)
{
  CVUserSuperLUMem cv_user_lmem;
  CVodeMem cvode_mem = (CVodeMem)cv_mem;

  if(cv_mem == NULL) {
    // TODO: issue cvode memory error
    printf("CVode Memory pointer is NULL\n"); fflush(stdout);
    return -1;
  }
  cv_user_lmem = (CVUserSuperLUMem) cvode_mem->cv_lmem;
  (*num_jac) = cv_user_lmem->num_jac;
  return CV_SUCCESS;
}

int CVUserSuperLUGetNumJacFactors(void *cv_mem,
                                 int *num_factor)
{
  CVUserSuperLUMem cv_user_lmem;
  CVodeMem cvode_mem = (CVodeMem)cv_mem;

  if(cv_mem == NULL) {
    // TODO: issue cvode memory error
    printf("CVode Memory pointer is NULL\n"); fflush(stdout);
    return -1;
  }
  cv_user_lmem = (CVUserSuperLUMem) cvode_mem->cv_lmem;
  (*num_factor) = cv_user_lmem->num_factor;
  return CV_SUCCESS;
}
int CVUserSuperLUGetNumJacSolves(void *cv_mem,
                                 int *num_solve)
{
  CVUserSuperLUMem cv_user_lmem;
  CVodeMem cvode_mem = (CVodeMem)cv_mem;

  if(cv_mem == NULL) {
    // TODO: issue cvode memory error
    printf("CVode Memory pointer is NULL\n"); fflush(stdout);
    return -1;
  }
  cv_user_lmem = (CVUserSuperLUMem) cvode_mem->cv_lmem;
  (*num_solve) = cv_user_lmem->num_solve;
  return CV_SUCCESS;
}
