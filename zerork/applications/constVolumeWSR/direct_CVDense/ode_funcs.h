#ifndef ODE_FUNCS_H
#define ODE_FUNCS_H

#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros

#include <zerork/mechanism.h>

typedef struct {
  double function_time;
  double jacobian_time;
  double lu_factor_time;
  double lu_solve_time;
} ODETimerStruct;

typedef struct {
  zerork::mechanism *mechp;
  double *rate_workspace;
  double concentration_ref;
  double inv_concentration_ref;
  double temperature_ref;
  double inv_temperature_ref;
  int num_roots;
  double *temperature_roots;
  ODETimerStruct timer;
} ODEParams_CMT;


ODEParams_CMT *AllocateODEParams_CMT(zerork::mechanism *mech,
                                     const double conc_ref,
                                     const double temp_ref,
                                     const int n_roots,
                                     const double temp_roots[]);
void ResetODEParams_CMT(const double conc_ref,
                        ODEParams_CMT *w);
void FreeODEParams_CMT(ODEParams_CMT *w);

int ConstantVolumeWSR_CMT(realtype t, N_Vector y, N_Vector ydot,
			  void *user_data);

int TempRootFunc_CMT(realtype t, N_Vector y, realtype *rootFunc,
		     void *user_data);
void ResetODETimerStruct(ODETimerStruct *w);
void AddODETimerStruct(ODETimerStruct *dest, ODETimerStruct *src);
#endif
