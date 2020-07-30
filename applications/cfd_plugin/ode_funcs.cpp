#include "ode_funcs.h"
#include "cv_param_sparse.h"
#include "utility_funcs.h"

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

  // set concentration via density and mass fraction
  cvp->mech->getCfromVY(cvp->invDens,y_ptr,cvp->conc);
  // set temperature
  Temp = y_ptr[cvp->nSpc]*cvp->Tref;
  
  // compute the molar production rates at the current state (aka wdot)
  cvp->mech->getReactionRates(Temp,cvp->conc,cvp->netProd,cvp->createRate,
			      cvp->destroyRate,cvp->fwdROP);

  cvp->mech->getIntEnergy_RT(Temp,cvp->Energy);
  cvp->meanCvMass=cvp->mech->getMassCvFromTY(Temp,y_ptr,cvp->CvMass);

  // ydot = [kmol/m^3/s] * [kg/kmol] * [m^3/kg] = [(kg spec j)/(kg mix)/s]
  for(j=0; j<cvp->nSpc; j++)
    {ydot_ptr[j]=(cvp->netProd[j])*(cvp->molWt[j])*(cvp->invDens);}

  Esum=0.0;
  for(j=0; j<cvp->nSpc; j++)
    {Esum+=(cvp->Energy[j])*(cvp->netProd[j]);}
 
  Esum*=cvp->mech->getGasConstant()*y_ptr[cvp->nSpc]*(cvp->invDens)
    /(cvp->meanCvMass);
  cvp->dTemp_dt=-Esum;

  ydot_ptr[cvp->nSpc]=cvp->dTemp_dt;

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

