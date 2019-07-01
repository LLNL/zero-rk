#include "ode_funcs.h"
#include "cv_param_sparse.h"
#include "utility_funcs.h"


#ifdef ZERORK_NEG_CONC_CHECK
#define ATOL_SAFETY_FACTOR 1.0e6
#endif

int const_vol_wsr(realtype t, N_Vector y, N_Vector ydot,
			 void *user_data)
{
  cv_param *cvp=(cv_param *)user_data;
  int nState = cvp->nSpc + 1;
  double * y_ptr = NV_DATA_S(y);
  double * ydot_ptr = NV_DATA_S(ydot);
  double Esum,Temp;
  int j;
  double startTime=getHighResolutionTime();

#ifdef ZERORK_NEG_CONC_CHECK
  //check for growing negative species
  // and try to get CVODE to recover
  const double check_val = -cvp->atol*ATOL_SAFETY_FACTOR;
  for(j=0; j<cvp->nSpc; ++j)
  {
    if(y_ptr[j] < check_val)
    {
      printf("WARNING: Significantly negative species in "
             "Zero-RK RHS.  Consider reducing tolerances.\n");
      return 1;
    }
  }
#endif

  // set concentration via density and mass fraction
  cvp->mech->getCfromVY(cvp->invDens[cvp->currReactor],y_ptr,cvp->conc);
  // set temperature
  Temp = y_ptr[cvp->nSpc]*cvp->Tref;
  
  // compute the molar production rates at the current state (aka wdot)
  cvp->mech->getReactionRates(Temp,cvp->conc,cvp->netProd,cvp->createRate,
			      cvp->destroyRate,cvp->fwdROP);

  cvp->mech->getIntEnergy_RT(Temp,cvp->Energy);
  cvp->meanCvMass[cvp->currReactor]=cvp->mech->getMassCvFromTY(Temp,y_ptr,cvp->CvMass);
  
  // ydot = [kmol/m^3/s] * [kg/kmol] * [m^3/kg] = [(kg spec j)/(kg mix)/s]
  for(j=0; j<cvp->nSpc; j++)
//    {NV_Ith_S(ydot,j)=(cvp->netProd[j])*(cvp->molWt[j])*(cvp->invDens[cvp->currReactor]);}
    {ydot_ptr[j]=(cvp->netProd[j])*(cvp->molWt[j])*(cvp->invDens[cvp->currReactor]);}

  Esum=0.0;
  for(j=0; j<cvp->nSpc; j++)
    {Esum+=(cvp->Energy[j])*(cvp->netProd[j]);}
 
  Esum*=cvp->mech->getGasConstant()*y_ptr[cvp->nSpc]*(cvp->invDens[cvp->currReactor])
    /(cvp->meanCvMass[cvp->currReactor]);
  cvp->dTemp_dt[cvp->currReactor]=-Esum;

  ydot_ptr[cvp->nSpc]=cvp->dTemp_dt[cvp->currReactor];

  (cvp->nFunc)++;
  cvp->funcTime += getHighResolutionTime() - startTime;
  return 0;
}


int const_dpdt_wsr(realtype t, N_Vector y, N_Vector ydot,
			 void *user_data)
{
  cv_param *cvp=(cv_param *)user_data;
  int nState = cvp->nSpc + 1;
  double * y_ptr = NV_DATA_S(y);
  double * ydot_ptr = NV_DATA_S(ydot);
  double Esum,Temp,invDensCurr,currPress;
  int j;
  double startTime=getHighResolutionTime();
  //N.B. Energy is Enthalpy and Cv is Cp in this routine.


#ifdef ZERORK_NEG_CONC_CHECK
  //check for growing negative species
  // and try to get CVODE to recover
  const double check_val = -cvp->atol*ATOL_SAFETY_FACTOR;
  for(j=0; j<cvp->nSpc; ++j)
  {
    if(y_ptr[j] < check_val)
    {
      printf("WARNING: Significantly negative species in "
             "Zero-RK RHS.  Consider reducing tolerances.\n");
      return 1;
    }
  }
#endif


  //ASSUME : t0 == 0.0;
#ifdef CONVERGE_CONST_PRESS
  currPress = cvp->Press[cvp->currReactor];
#else
#ifdef CONVERGE_DPDT_UNITS
  currPress = cvp->Press[cvp->currReactor]+(t)*cvp->dpdt[cvp->currReactor]/10.0;
#else
  currPress = cvp->Press[cvp->currReactor]+(t)*cvp->dpdt[cvp->currReactor];
#endif
#endif

  // set temperature
  Temp = y_ptr[cvp->nSpc]*cvp->Tref;

  invDensCurr = cvp->mech->getDensityFromTPY(Temp,currPress,y_ptr);
  invDensCurr = 1.0/invDensCurr;

  // set concentration via density and mass fraction
  cvp->mech->getCfromVY(invDensCurr,y_ptr,cvp->conc);

  // compute the molar production rates at the current state (aka wdot)
  cvp->mech->getReactionRates(Temp,cvp->conc,cvp->netProd,cvp->createRate,
			      cvp->destroyRate,cvp->fwdROP);

  cvp->mech->getEnthalpy_RT(Temp,cvp->Energy);
  cvp->meanCvMass[cvp->currReactor]=cvp->mech->getMassCpFromTY(Temp,y_ptr,cvp->CvMass);

  // ydot = [kmol/m^3/s] * [kg/kmol] * [m^3/kg] = [(kg spec j)/(kg mix)/s]
  for(j=0; j<cvp->nSpc; j++)
    {ydot_ptr[j]=(cvp->netProd[j])*(cvp->molWt[j])*(invDensCurr);}

  Esum=0.0;
  for(j=0; j<cvp->nSpc; j++)
    {Esum+=(cvp->Energy[j])*(cvp->netProd[j]);}

  Esum*=cvp->mech->getGasConstant()*y_ptr[cvp->nSpc]*(invDensCurr)
    /(cvp->meanCvMass[cvp->currReactor]);
  //CONVERGE dpdt term
#ifdef CONVERGE_DPDT_UNITS
  Esum-=(cvp->dpdt[cvp->currReactor]*invDensCurr/10.0
           /(cvp->meanCvMass[cvp->currReactor]*cvp->Tref));
#else
  Esum-=(cvp->dpdt[cvp->currReactor]*invDensCurr
           /(cvp->meanCvMass[cvp->currReactor]*cvp->Tref));
#endif
  //Sign change
  cvp->dTemp_dt[cvp->currReactor]=-Esum;

  ydot_ptr[cvp->nSpc]=cvp->dTemp_dt[cvp->currReactor];

  (cvp->nFunc)++;
  cvp->funcTime += getHighResolutionTime() - startTime;
  return 0;
}


int tempRootFunc(realtype t, N_Vector y, realtype *rootFunc,
                void *user_data)
{
  cv_param *cvp=(cv_param *)user_data;
  int nSize = cvp->sparseMtx->nSize; //single reactor system size (nSpc + 1)
  int nSpc = nSize-1;

  rootFunc[0]=(cvp->Tinit+cvp->deltaTign)-(NV_Ith_S(y,cvp->nSpc)*cvp->Tref);
  return 0;
}

