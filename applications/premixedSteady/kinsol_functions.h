#ifndef KINSOL_FUNCTIONS_H_
#define KINSOL_FUNCTIONS_H_

#include <kinsol/kinsol.h>            // prototypes for KINSOL fcts. and consts.
#ifdef ZERORK_MPI
#include <nvector/nvector_parallel.h> // serial N_Vector types, fcts., and macros
#else
#include <nvector/nvector_serial.h>
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

void ErrorFunction(int error_code,
                   const char *module,
                   const char *function,
                   char *msg,
                   void *user_data);

// Approximate factorization preconditioner solve
int AFSolve(double solution[],
            void *user_data);


typedef struct
{
  int rxnId;
  double relSens;
  double relSensAbs;
} rxnSens_t;

int compare_rxnSens_t(const void *A, const void *B);

// Flame speed reaction sensitivity analysis
int SensitivityAnalysis(N_Vector y,
                        N_Vector ydot,
			void *user_data);

#endif
