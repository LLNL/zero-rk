#ifndef CVODE_FUNCTIONS_H_
#define CVODE_FUNCTIONS_H_

#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_parallel.h> // serial N_Vector types, fcts., and macros

int ConstPressureFlame(realtype t,
		       N_Vector y,
		       N_Vector ydot,
		       void *user_data);


int ConstPressureFlameComm(realtype t,
			   int nlocal,
			   N_Vector y,
			   void *user_data);

int ConstPressureFlameLocal(int nlocal,
			    realtype t,
			    N_Vector y,
			    N_Vector ydot,
			    void *user_data);

#if defined SUNDIALS2
int ReactorPreconditionerSetup(realtype t,// [in] ODE system time
                               N_Vector y,      // [in] ODE state vector
                               N_Vector ydot,   // [in] ODE state derivative
                               booleantype jok,
			       booleantype *new_j,
			       realtype gamma,
		               void *params,
                               N_Vector tmp1,
                               N_Vector tmp2,
                               N_Vector tmp3);

int ReactorPreconditionerSolve(realtype t,      // [in] ODE system time
                               N_Vector y,      // [in] ODE state vector
                               N_Vector ydot,   // [in] ODE state derivative
                               N_Vector r,      // [in] jacobian rhs
                               N_Vector z,      // [out]
			       realtype gamma,
			       realtype delta,
			       int lr,
		               void *params,
                               N_Vector tmp);

int ReactorBBDSetup(realtype t,// [in] ODE system time
                    N_Vector y,      // [in] ODE state vector
                    N_Vector ydot,   // [in] ODE state derivative
                    booleantype jok,
                    booleantype *new_j,
                    realtype gamma,
                    void *params,
                    N_Vector tmp1,
                    N_Vector tmp2,
                    N_Vector tmp3);

int ReactorBBDSolve(realtype t,      // [in] ODE system time
                    N_Vector y,      // [in] ODE state vector
                    N_Vector ydot,   // [in] ODE state derivative
                    N_Vector r,      // [in] jacobian rhs
                    N_Vector z,      // [out]
                    realtype gamma,
                    realtype delta,
                    int lr,
                    void *params,
                    N_Vector tmp);
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorPreconditionerSetup(realtype t,// [in] ODE system time
                               N_Vector y,      // [in] ODE state vector
                               N_Vector ydot,   // [in] ODE state derivative
                               booleantype jok,
			       booleantype *new_j,
			       realtype gamma,
		               void *params);

int ReactorPreconditionerSolve(realtype t,      // [in] ODE system time
                               N_Vector y,      // [in] ODE state vector
                               N_Vector ydot,   // [in] ODE state derivative
                               N_Vector r,      // [in] jacobian rhs
                               N_Vector z,      // [out]
			       realtype gamma,
			       realtype delta,
			       int lr,
		               void *params);

int ReactorBBDSetup(realtype t,// [in] ODE system time
                    N_Vector y,      // [in] ODE state vector
                    N_Vector ydot,   // [in] ODE state derivative
                    booleantype jok,
                    booleantype *new_j,
                    realtype gamma,
                    void *params);

int ReactorBBDSolve(realtype t,      // [in] ODE system time
                    N_Vector y,      // [in] ODE state vector
                    N_Vector ydot,   // [in] ODE state derivative
                    N_Vector r,      // [in] jacobian rhs
                    N_Vector z,      // [out]
                    realtype gamma,
                    realtype delta,
                    int lr,
                    void *params);
#endif

#endif
