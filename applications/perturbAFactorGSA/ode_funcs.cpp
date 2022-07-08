#include "ode_funcs.h"
#include "cv_param_sparse.h"
#include "utility_funcs.h"

#include "utilities.h"

using zerork::getHighResolutionTime;

// int const_vol_wsr(realtype t, N_Vector y, N_Vector ydot,
// 			 void *user_data)
// {
//   cv_param *cvp=(cv_param *)user_data;
//   double *mfp=NV_DATA_S(y); // caution: assumes realtype == double
//   double Esum,Temp;
//   int j;
//   double startTime=getHighResolutionTime();
  
//   // set concentration via density and mass fraction
//   cvp->mechPtr->getCfromVY(cvp->invDens,mfp,cvp->conc);
//   // set temperature
//   Temp = NV_Ith_S(y,cvp->nSpc)*cvp->Tref;
  
//   // compute the molar production rates at the current state (aka wdot)
//   cvp->mechPtr->getReactionRates(Temp,cvp->conc,cvp->netProd,cvp->createRate,
// 			      cvp->destroyRate,cvp->fwdROP);

//   cvp->mechPtr->getIntEnergy_RT(Temp,cvp->Energy);
//   cvp->meanCvMass=cvp->mechPtr->getMassCvFromTY(Temp,mfp,cvp->CvMass);
  
//   // ydot = [kmol/m^3/s] * [kg/kmol] * [m^3/kg] = [(kg spec j)/(kg mix)/s]
//   for(j=0; j<cvp->nSpc; j++)
//     {NV_Ith_S(ydot,j)=(cvp->netProd[j])*(cvp->molWt[j])*(cvp->invDens);}

//   Esum=0.0;
//   for(j=0; j<cvp->nSpc; j++)
//     {Esum+=(cvp->Energy[j])*(cvp->netProd[j]);}
 
//   Esum*=cvp->mechPtr->getGasConstant()*NV_Ith_S(y,cvp->nSpc)*(cvp->invDens)
//     /(cvp->meanCvMass);
//   cvp->dTemp_dt=-Esum;

//   NV_Ith_S(ydot,cvp->nSpc)=cvp->dTemp_dt;

//   (cvp->nFunc)++;
//   cvp->funcTime += getHighResolutionTime() - startTime;
//   return 0;
// }

int const_vol_wsr_perturb(realtype t, N_Vector y, N_Vector ydot,
			  void *user_data)
{
  cv_param *cvp=(cv_param *)user_data;
  double *mfp=NV_DATA_S(y); // caution: assumes realtype == double
  double Esum,Temp;
  int j;
  double startTime=getHighResolutionTime();
  
  // set concentration via density and mass fraction
  cvp->mechPtr->getCfromVY(cvp->invDens,mfp,cvp->conc);
  // set temperature
  Temp = NV_Ith_S(y,cvp->nSpc)*cvp->Tref;
  
  // compute the molar production rates at the current state (aka wdot)
  cvp->mechPtr->getReactionRates_perturbROP(Temp,
                                            cvp->conc,
                                            cvp->ropMultiplierPtr,
                                            cvp->netProd,
                                            cvp->createRate,
			                    cvp->destroyRate,
                                            cvp->fwdROP);

  cvp->mechPtr->getIntEnergy_RT(Temp,cvp->Energy);
  cvp->meanCvMass=cvp->mechPtr->getMassCvFromTY(Temp,mfp,cvp->CvMass);
  
  // ydot = [kmol/m^3/s] * [kg/kmol] * [m^3/kg] = [(kg spec j)/(kg mix)/s]
  for(j=0; j<cvp->nSpc; j++)
    {NV_Ith_S(ydot,j)=(cvp->netProd[j])*(cvp->molWt[j])*(cvp->invDens);}

  Esum=0.0;
  for(j=0; j<cvp->nSpc; j++)
    {Esum+=(cvp->Energy[j])*(cvp->netProd[j]);}
 
  Esum*=cvp->mechPtr->getGasConstant()*NV_Ith_S(y,cvp->nSpc)*(cvp->invDens)
    /(cvp->meanCvMass);
  cvp->dTemp_dt=-Esum;

  NV_Ith_S(ydot,cvp->nSpc)=cvp->dTemp_dt;

  (cvp->nFunc)++;
  cvp->funcTime += getHighResolutionTime() - startTime;
  return 0;
}

int tempRootFunc(realtype t, N_Vector y, realtype *rootFunc,
		 void *user_data)
{
  int j;
  cv_param *cvp=(cv_param *)user_data;
  
  for(j=0; j<cvp->nIdtTemp; j++) {
    rootFunc[j]=cvp->redTempRoot[j]-NV_Ith_S(y,cvp->nSpc);
  }
  return 0;
}

