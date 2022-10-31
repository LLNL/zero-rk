#include "ode_funcs.h"
#include "cv_param_sparse.h"
#include "utility_funcs.h"

using zerork::getHighResolutionTime;

int const_vol_wsr(realtype t, N_Vector y, N_Vector ydot,
			 void *user_data)
{
  cv_param *cvp=(cv_param *)user_data;
  double *mfp=NV_DATA_S(y); // caution: assumes realtype == double
  double Esum,Temp;
  int j;
  double startTime=getHighResolutionTime();

  // set concentration via density and mass fraction
  cvp->mech->getCfromVY(cvp->invDens,mfp,cvp->conc);

  // set temperature
  Temp = NV_Ith_S(y,cvp->nSpc)*cvp->Tref;

  // compute the molar production rates at the current state (aka wdot)
  cvp->mech->getReactionRates(Temp,cvp->conc,cvp->netProd,cvp->createRate,
			      cvp->destroyRate,cvp->fwdROP);

  cvp->mech->getIntEnergy_RT(Temp,cvp->Energy);
  cvp->meanCvMass=cvp->mech->getMassCvFromTY(Temp,mfp,cvp->CvMass);

  // ydot = [kmol/m^3/s] * [kg/kmol] * [m^3/kg] = [(kg spec j)/(kg mix)/s]
  for(j=0; j<cvp->nSpc; j++)
    {NV_Ith_S(ydot,j)=(cvp->netProd[j])*(cvp->molWt[j])*(cvp->invDens);}

  Esum=0.0;
  for(j=0; j<cvp->nSpc; j++)
    {Esum+=(cvp->Energy[j])*(cvp->netProd[j]);}

  Esum*=cvp->mech->getGasConstant()*NV_Ith_S(y,cvp->nSpc)*(cvp->invDens)
    /(cvp->meanCvMass);
  cvp->dTemp_dt=-Esum;

  NV_Ith_S(ydot,cvp->nSpc)=cvp->dTemp_dt;

  (cvp->nFunc)++;
  cvp->funcTime += getHighResolutionTime() - startTime;
  return 0;
}

int const_vol_wsr_Sij(realtype t, N_Vector y, N_Vector ydot, double Sij[],
			 void *user_data)
{
  cv_param *cvp=(cv_param *)user_data;
  double *mfp=NV_DATA_S(y); // caution: assumes realtype == double
  double Esum,Temp;
  int j;
  double startTime=getHighResolutionTime();

  // set concentration via density and mass fraction
  cvp->mech->getCfromVY(cvp->invDens,mfp,cvp->conc);

  // set temperature
  Temp = NV_Ith_S(y,cvp->nSpc)*cvp->Tref;

  // compute the molar production rates at the current state (aka wdot)
  cvp->mech->getReactionRates(Temp,cvp->conc,cvp->netProd,cvp->createRate,
			      cvp->destroyRate,cvp->fwdROP);

  cvp->mech->getIntEnergy_RT(Temp,cvp->Energy);
  cvp->meanCvMass=cvp->mech->getMassCvFromTY(Temp,mfp,cvp->CvMass);

  // ydot = [kmol/m^3/s] * [kg/kmol] * [m^3/kg] = [(kg spec j)/(kg mix)/s]
  for(j=0; j<cvp->nSpc; j++)
    {NV_Ith_S(ydot,j)=(cvp->netProd[j])*(cvp->molWt[j])*(cvp->invDens);}

  Esum=0.0;
  for(j=0; j<cvp->nSpc; j++)
    {Esum+=(cvp->Energy[j])*(cvp->netProd[j]);}

  Esum*=cvp->mech->getGasConstant()*NV_Ith_S(y,cvp->nSpc)*(cvp->invDens)
    /(cvp->meanCvMass);
  cvp->dTemp_dt=-Esum;

  NV_Ith_S(ydot,cvp->nSpc)=cvp->dTemp_dt;

  int nStep = cvp->mech->getNumSteps();
  int nRxn = cvp->mech->getNumReactions();
  int i,k;
  int nreac, nprod, fwdId, revId;
  std::vector<double> w_ij;
  w_ij.assign(cvp->nSpc*nRxn, 0.0);

  // Loop over reactions to compute w_ij
  for(k=0; k<nRxn; k++) {

    // Get forward and reverse step index
    fwdId = cvp->mech->getStepIdxOfRxn(k,1);
    revId = cvp->mech->getStepIdxOfRxn(k,-1);

    // Forward reaction
    nreac = cvp->mech->getOrderOfStep(fwdId);
    nprod = cvp->mech->getNumProductsOfStep(fwdId);

    for(i=0; i<nreac; i++) {
      j = cvp->mech->getSpecIdxOfStepReactant(fwdId,i);
      w_ij[k*cvp->nSpc + j] -= cvp->fwdROP[fwdId];
    }

    for(i=0; i<nprod; i++) {
      j = cvp->mech->getSpecIdxOfStepProduct(fwdId,i);
      w_ij[k*cvp->nSpc + j] += cvp->fwdROP[fwdId];
    }

    // Reverse reaction
    if(revId >= 0 && revId < nStep) {
      nreac = cvp->mech->getOrderOfStep(revId);
      nprod = cvp->mech->getNumProductsOfStep(revId);

      for(i=0; i<nreac; i++) {
        j = cvp->mech->getSpecIdxOfStepReactant(revId,i);
        w_ij[k*cvp->nSpc + j] -= cvp->fwdROP[revId];
      }

      for(i=0; i<nprod; i++) {
        j = cvp->mech->getSpecIdxOfStepProduct(revId,i);
        w_ij[k*cvp->nSpc + j] += cvp->fwdROP[revId];
      }
    }
  }

  // Loop over reactions to compute S_ij
  for(k=0; k<nRxn; k++) {
    // Species
    for(j=0; j<cvp->nSpc; j++) {
      Sij[k*(cvp->nSpc+1)+j] = w_ij[k*cvp->nSpc+j]*
        (cvp->molWt[j])*(cvp->invDens);
    }

    // Temperature
    Esum=0.0;
    for(j=0; j<cvp->nSpc; j++) {
      Esum+=(cvp->Energy[j])*(w_ij[k*cvp->nSpc+j]);
    }
    Esum*=cvp->mech->getGasConstant()*mfp[cvp->nSpc]*(cvp->invDens)/
      (cvp->meanCvMass);
    Sij[k*(cvp->nSpc+1) + cvp->nSpc] = -Esum;
  }

  (cvp->nFunc)++;
  cvp->funcTime += getHighResolutionTime() - startTime;
  return 0;
}

int tempRootFunc(realtype t, N_Vector y, realtype *rootFunc,
		 void *user_data)
{
  cv_param *cvp=(cv_param *)user_data;

  for(int i = 0; i < cvp->numTempRoots; ++i) {
      rootFunc[i]=(cvp->tempRoots[i])-(NV_Ith_S(y,cvp->nSpc)*cvp->Tref);
  }
  return 0;
}
