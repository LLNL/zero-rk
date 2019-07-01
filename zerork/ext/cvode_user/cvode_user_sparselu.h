#ifndef CVODE_USER_SUPERLU_H_
#define CVODE_USER_SUPERLU_H_

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros
#include <cvode/cvode_impl.h>
#include <slu_ddefs.h>

  typedef struct CVUserSuperLUMemRec{

    int num_equations;
    int num_nonzeros;
 
    int *row_index;
    int *column_sum;
    double *J_user;
    double *M_user;
    double *b_user;
    double *x_user;
    double *sparse_identity;
    int (*UserJacobianFunction)(int N,
                                realtype t,
                                N_Vector y,
                                N_Vector fy,
                                double user_jac[],
                                void *user_data,
                                N_Vector tmp1,
                                N_Vector tmp2,
                                N_Vector tmp3);
    

    // SuperLU data
    int permutationType;
    char equed[1];
    SuperMatrix Mslu;
    SuperMatrix Lslu;
    SuperMatrix Uslu;
    SuperMatrix Bslu;
    SuperMatrix Xslu;
    int *row_permutation;
    int *col_permutation;
    int *col_elim_tree; // etree
    double *Rvec_internal;
    double *Cvec_internal;
    superlu_options_t optionSLU;
    SuperLUStat_t statSLU;
    mem_usage_t  mem_usage;

    // solver counter, flags
    int num_jac;
    int num_factor;
    int num_solve;
    int num_steps_last_jac;
    int last_flag;

  } *CVUserSuperLUMem;

  // Four basic functions to provide an alternate linear solver
  //  1. Initialization
  int LInitUserSuperLU(CVodeMem cv_mem);

  //  2. Setup - prepares the linear solver for subsequent calls to LSolve 
  int LSetupUserSuperLU(CVodeMem cv_mem,
                        int convfail,
                        N_Vector ypred,
                        N_Vector fpred,
                        booleantype *jcurPtr,
                        N_Vector vtemp1,
                        N_Vector vtemp2,
                        N_Vector vtemp3);

  // 3. Solve - solve the linear system
  int LSolveUserSuperLU(CVodeMem cv_mem,
                        N_Vector b,
                        N_Vector weight,
                        N_Vector ycur,
                        N_Vector fcur);

  // 4. Free any allocated memory
  void LFreeUserSuperLU(CVodeMem cv_mem);

  int CVUserSuperLU(void *cv_mem,
                    const int num_eq,
                    const int num_nz,
                    const int row_id[],
		    const int col_start[],
                    int (*jac_func_ptr)(int,realtype, N_Vector, N_Vector,
                                        double *, void *, N_Vector, N_Vector,
                                        N_Vector));

  int HasCompleteDiagonal(const int num_eq,
                          const int num_nz,
                          const int row_id[],
		          const int col_start[],
                          int diag_sparse_id[]);

  int CVUserSuperLUGetNumJacEvals(void *cv_mem,
                                  int *num_jac);
  int CVUserSuperLUGetNumJacFactors(void *cv_mem,
                                    int *num_factor);
  int CVUserSuperLUGetNumJacSolves(void *cv_mem,
                                   int *num_solve);

  

#ifdef __cplusplus
}
#endif

#endif
