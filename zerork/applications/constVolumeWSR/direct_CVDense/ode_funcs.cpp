#include <stdlib.h>
#include <stdio.h>


#include <derivative/const_vol.h>

#include "utility_funcs.h"

#include "ode_funcs.h"


int ConstantVolumeWSR_CMT(realtype t, N_Vector y, N_Vector ydot,
			  void *user_data)
{
  int j, return_flag;
  double start_time=getHighResolutionTime();
  ODEParams_CMT *params    = (ODEParams_CMT *)user_data;
  const int num_species    = params->mechp->getNumSpecies();
  double temperature;
  double mixture_conc_rate = 0.0;
  double temperature_rate;
  double *concentration = &params->rate_workspace[0];
  double *net_rate = NV_DATA_S(ydot);
  double *func_work     = &params->rate_workspace[num_species];

  for(j=0; j<num_species; ++j) {
    concentration[j] = NV_Ith_S(y,j)*params->concentration_ref;
  }
  temperature = NV_Ith_S(y,num_species+1)*params->temperature_ref;

  return_flag = zerork::ConstVolume_CT(params->mechp,
                                    concentration,
                                    temperature,
				    func_work,
                                    net_rate,
                                    &temperature_rate);    
  for(j=0; j<num_species; ++j) {
    net_rate[j]*=params->inv_concentration_ref;                             
    mixture_conc_rate += net_rate[j];
  }
  NV_Ith_S(ydot,num_species)   = mixture_conc_rate;
  NV_Ith_S(ydot,num_species+1) = temperature_rate*params->inv_temperature_ref; 
  params->timer.function_time+=getHighResolutionTime()-start_time;
  return return_flag;
}

int TempRootFunc_CMT(realtype t, N_Vector y, realtype *rootFunc,
		     void *user_data)
{
  int j;
  ODEParams_CMT *cv_params=(ODEParams_CMT *)user_data;
  const int num_roots = cv_params->num_roots;
  const int num_species = cv_params->mechp->getNumSpecies();

  for(j=0; j<num_roots; ++j) {
    rootFunc[j]=(cv_params->temperature_roots[j])-
      (NV_Ith_S(y,num_species+1)*cv_params->temperature_ref);
  }
  return 0;
}

ODEParams_CMT *AllocateODEParams_CMT(zerork::mechanism *mech,
                                     const double conc_ref,
                                     const double temp_ref,
                                     const int num_temp_roots,
                                     const double temp_roots[])
{
  int j;
  const int num_species = mech->getNumSpecies();
  ODEParams_CMT *w;
  
  w =(ODEParams_CMT *)malloc(sizeof(ODEParams_CMT));
  // TODO: add memory allocation checks

  // basic assignment
  w->mechp = mech;
  w->concentration_ref = conc_ref;
  w->inv_concentration_ref = 1.0/conc_ref;
  w->temperature_ref = temp_ref;
  w->inv_temperature_ref = 1.0/temp_ref;
  w->num_roots = num_temp_roots;

  printf("Number of species: %d\n", w->mechp->getNumSpecies());
  // initialization of timer
  ResetODETimerStruct(&w->timer);  
 
  // allocate arrays
  w->rate_workspace    = (double *)malloc(sizeof(double)*2*num_species);
  w->temperature_roots = (double *)malloc(sizeof(double)*num_temp_roots);

  // assign temperature roots
  for(j=0; j<num_temp_roots; ++j) {
    w->temperature_roots[j] = temp_roots[j];
  }

  return w;
}
void ResetODEParams_CMT(const double conc_ref,
                        ODEParams_CMT *w)
{
  w->concentration_ref = conc_ref;
  w->inv_concentration_ref = 1.0/conc_ref;
  ResetODETimerStruct(&w->timer);
}

void FreeODEParams_CMT(ODEParams_CMT *w)
{
  free(w->temperature_roots);
  free(w->rate_workspace);
  free(w);
}

void ResetODETimerStruct(ODETimerStruct *w)
{
  w->function_time  = 0.0;
  w->jacobian_time  = 0.0;
  w->lu_factor_time = 0.0;
  w->lu_solve_time  = 0.0;
}
void AddODETimerStruct(ODETimerStruct *dest, ODETimerStruct *src)
{
  dest->function_time  += src->function_time;
  dest->jacobian_time  += src->jacobian_time;
  dest->lu_factor_time += src->lu_factor_time;
  dest->lu_solve_time  += src->lu_solve_time;
}
