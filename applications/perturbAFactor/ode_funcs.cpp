#include "ode_funcs.h"
#include "cv_param_sparse.h"
#include "utility_funcs.h"

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
  const int num_species = cvp->nSpc;
  double *mfp=NV_DATA_S(y); // caution: assumes realtype == double
  double *derivative = NV_DATA_S(ydot);
  double Esum,Temp;
  int j;
  double startTime=getHighResolutionTime();
  
  // set concentration via density and mass fraction
  cvp->mechPtr->getCfromVY(cvp->invDens,mfp,cvp->conc);
  // set temperature
  Temp = NV_Ith_S(y,num_species)*cvp->Tref;
  
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
  for(j=0; j<num_species; j++) {
    derivative[j]=(cvp->netProd[j])*(cvp->molWt[j])*(cvp->invDens);
  }
  Esum=0.0;
  for(j=0; j<num_species; j++) {
    Esum+=(cvp->Energy[j])*(cvp->netProd[j]);
  }
 
  Esum*=cvp->mechPtr->getGasConstant()*NV_Ith_S(y,num_species)*(cvp->invDens)
    /(cvp->meanCvMass);
  cvp->dTemp_dt = -Esum;

  derivative[num_species]=cvp->dTemp_dt;

  (cvp->nFunc)++;
  cvp->funcTime += getHighResolutionTime() - startTime;
  return 0;
}


// tempRootFunc[0] = (T_idt[0] - T)/T_ref,  want zero-crossing from + to -
//             [:]
// tempRootFunc[num_idt_temperatures_] = ydot[track_species_max_id_[0]
//                                       want zero-crossing from + to -
//                                       for maxima
int tempRootFunc(realtype t, 
                 N_Vector y, 
                 realtype *rootFunc,
		 void *user_data)
{
  int j;
  int flag;
  cv_param *cvp=(cv_param *)user_data;
  long int num_cvode_steps;
  const int num_idt_temp = cvp->nIdtTemp;
  const int num_track_species_max =  cvp->num_track_species_max_;
  double *derivative = NV_DATA_S(cvp->state_workspace_);
  flag = CVodeGetNumSteps(cvp->cvodeMemPtr, &num_cvode_steps);
  if(flag < 0) {
    return flag;
  }

  for(j=0; j<num_idt_temp; j++) {
    rootFunc[j]=cvp->redTempRoot[j]-NV_Ith_S(y,cvp->nSpc);
  }

  if(0) {
    // try to compute the rate of change of the mass fraction for the
    // root function to determine the maximum value
    if(num_track_species_max > 0) {

      if(num_cvode_steps > 0) {
        // get the value of dy/dt at time t
        flag = CVodeGetDky(cvp->cvodeMemPtr,
                           t,
                           1, // derivative order
                          cvp->state_workspace_);
      } else {
        // compute directly
        flag = const_vol_wsr_perturb(t,y,cvp->state_workspace_,user_data);
      }
      if(flag < 0) {
        return flag;
      }
    }

    for(j=0; j<num_track_species_max; ++j) {
      int track_id = cvp->track_species_max_id_[j];
      rootFunc[j+num_idt_temp] = derivative[track_id];
    }
  } else {
    // set the root function to one for the tracked species so a root is 
    // never found
    for(j=0; j<num_track_species_max; ++j) {
      rootFunc[j+num_idt_temp] = 1.0;
    }
  }
  return 0;
}

// function has the same prototype as CVQuadRhsFn
// NV_Ith_S(yQdot, 0) = -\rho \sum u_{f,i}(T_{ref})\frac{dy_i}{dt} [J/m^3/s]
// NV_Ith_S(yQdot, 1) = -\rho \sum u_i(T)\frac{dy_i}{dt} [J/m^3/s]
int ChemicalHeatReleaseRate(realtype t, 
                            N_Vector y,
                            N_Vector yQdot,
                            void *user_data)
{
  cv_param *cvp=(cv_param *)user_data;
  double *state = NV_DATA_S(y);
  double *derivative = NV_DATA_S(cvp->state_workspace_);
  double *integrand = NV_DATA_S(yQdot);
  double *internal_energy = cvp->conc; // use the concentration workspace
  double *energy_of_formation = cvp->energy_of_formation_; // [J/kg/K]
  const int num_species = cvp->nSpc;
  const double temperature = state[num_species]*cvp->Tref;
  const double density = cvp->Dens;
  const double gas_constant = cvp->mechPtr->getGasConstant();
  int flag;

  // get the current internal energy of each species in [J/kg/K]
  cvp->mechPtr->getIntEnergy_RT(temperature, internal_energy);
  for(int j=0; j<num_species; ++j) {
    internal_energy[j] *= gas_constant*temperature*cvp->invMolWt[j]; 
  }
  // get the value of dy/dt at time t
  flag = CVodeGetDky(cvp->cvodeMemPtr,
                     t,
                     1, // derivative order
                     cvp->state_workspace_);
  //for(int j=0; j<num_species; ++j) {
  //  printf("dy[%d]/dt = %14.7e %14.7e (E_f = %14.7e, E = %14.7e)\n", 
  //         j, 
  //         NV_Ith_S(cvp->state_workspace_,j),
  //         derivative[j],
  //         energy_of_formation[j],
  //         internal_energy[j]);
  //}

  integrand[0] = 0.0;
  integrand[1] = 0.0;
  if(flag != CV_SUCCESS) {
    return flag;
  }
  for(int j=0; j<num_species; ++j) {
    integrand[0] += derivative[j]*energy_of_formation[j];
    integrand[1] += derivative[j]*internal_energy[j];
  }  
  integrand[0] *= -density;
  integrand[1] *= -density;
  //printf("# DEBUG:  %14.7e  %14.7e  %14.7e\n",t, integrand[0], integrand[1]);  

  return CV_SUCCESS;
}
