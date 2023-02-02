#ifndef KINSOL_FUNCTIONS_H_
#define KINSOL_FUNCTIONS_H_

#include <kinsol/kinsol.h>    // prototypes for KINSOL fcts. and consts.
#ifdef ZERORK_MPI
#include <nvector/nvector_parallel.h> // N_Vector type
#else
#include <nvector/nvector_serial.h> // N_Vector type
#endif

int ConstPressureFlame(N_Vector y,
		       N_Vector ydot,
		       void *user_data);

int ConstPressureFlameComm(int nlocal,
			   N_Vector y,
			   void *user_data);

int ConstPressureFlameLocal(int nlocal,
			    N_Vector y,
			    N_Vector ydot,
			    void *user_data);
#if defined SUNDIALS2
// SuperLU distributed
int ReactorBBDSetup(N_Vector y, // [in] ODE state vector
		    N_Vector yscale, // [in] ODE state scaler
		    N_Vector ydot, // [in] ODE state derivative
		    N_Vector ydotscale, // [in] ODE state derivative scaler
		    void *params,  // [in/out]
                    N_Vector tmp1, N_Vector tmp2);

int ReactorBBDSolve(N_Vector y, // [in] ODE state vector
		    N_Vector yscale, // [in] ODE state scaler
		    N_Vector ydot, // [in] ODE state derivative
		    N_Vector ydotscale, // [in] ODE state derivatie scaler
		    N_Vector vv, // [in] ??
		    void *params, // [in/out]
                    N_Vector tmp);
// SuperLU + ScaLapack Approximate Factorization
int ReactorAFSetup(N_Vector y, // [in] ODE state vector
		   N_Vector yscale, // [in] ODE state scaler
		   N_Vector ydot, // [in] ODE state derivative
		   N_Vector ydotscale, // [in] ODE state derivative scaler
		   void *params, // [in/out]
                   N_Vector tmp1, N_Vector tmp2);
int ReactorAFSolve(N_Vector y, // [in] ODE state vector
		   N_Vector yscale, // [in] ODE state scaler
		   N_Vector ydot, // [in] ODE state derivative
		   N_Vector ydotscale, // [in] ODE state derivatie scaler
		   N_Vector vv, // [in] ??
		   void *params, // [in/out]
                   N_Vector tmp);
#elif defined SUNDIALS3 || defined SUNDIALS4
// SuperLU distributed
int ReactorBBDSetup(N_Vector y, // [in] ODE state vector
		    N_Vector yscale, // [in] ODE state scaler
		    N_Vector ydot, // [in] ODE state derivative
		    N_Vector ydotscale, // [in] ODE state derivative scaler
		    void *params); // [in/out]

int ReactorBBDSolve(N_Vector y, // [in] ODE state vector
		    N_Vector yscale, // [in] ODE state scaler
		    N_Vector ydot, // [in] ODE state derivative
		    N_Vector ydotscale, // [in] ODE state derivatie scaler
		    N_Vector vv, // [in] ??
		    void *params); // [in/out]

// SuperLU + ScaLapack Approximate Factorization
int ReactorAFSetup(N_Vector y, // [in] ODE state vector
		   N_Vector yscale, // [in] ODE state scaler
		   N_Vector ydot, // [in] ODE state derivative
		   N_Vector ydotscale, // [in] ODE state derivative scaler
		   void *params); // [in/out]

int ReactorAFSolve(N_Vector y, // [in] ODE state vector
		   N_Vector yscale, // [in] ODE state scaler
		   N_Vector ydot, // [in] ODE state derivative
		   N_Vector ydotscale, // [in] ODE state derivatie scaler
		   N_Vector vv, // [in] ??
		   void *params); // [in/out]
#endif

// Approximate factorization preconditioner solve
int AFSolve(double solution[],
            void *user_data);

#endif
