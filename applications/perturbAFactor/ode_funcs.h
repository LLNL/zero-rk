#ifndef ODE_FUNCS_H
#define ODE_FUNCS_H

#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros

//int const_vol_wsr(realtype t, N_Vector y, N_Vector ydot,
//		  void *user_data);
int const_vol_wsr_perturb(realtype t, N_Vector y, N_Vector ydot,
		          void *user_data);

int tempRootFunc(realtype t, N_Vector y, realtype *rootFunc,
		 void *user_data);

// function has the same prototype as CVQuadRhsFn
// NV_Ith_S(yQdot, 0) = -\rho \sum u_{f,i}(T_{ref})\frac{dy_i}{dt} [J/m^3/s]
// NV_Ith_S(yQdot, 1) = -\rho \sum u_i(T)\frac{dy_i}{dt} [J/m^3/s]
int ChemicalHeatReleaseRate(realtype t, 
                            N_Vector y,
                            N_Vector yQdot,
                            void *user_data);

#endif
