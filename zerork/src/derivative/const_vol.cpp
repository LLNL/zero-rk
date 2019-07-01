#include "const_vol.h"

namespace zerork {

// ConstVolume_CT(...)
//
// Computes the time derivatives of the composition and thermodynamic state
// for a constant volume, homogeneous reactor based on the CT state
// specification (C = concentration, T = temperature).  Note that the sum
// of the concentrations relates to the second thermodynamic property needed
// to uniquely specify the state (either density or pressure).
//
// Inputs (unmodified):
//   zerork::mechanism *mech - mechanism object
//   const double conc[] - species concentration [kmol/m^3], 
//                         length=mech->numSpecies()
//   const double temp   - temperature [K]
//  
// Outputs:
//   double workspace[] - workspace array min length=mech->numSpecies()
//   double d_conc_dt[] - time derivative of the concentrations [kmol/m^3/s]
//                        length=mech->numSpecies()
//   double *d_temp_dt  - pointer to time derivative of temperature [K/s]
//
// Return:
//   int errorFlag      - returns 0 if no error, or -1 if there is an error
int ConstVolume_CT(zerork::mechanism *mech,
                   const double conc[],
                   const double temp,
                   double workspace[],
                   double d_conc_dt[],
                   double *d_temp_dt)
{
  int j;
  double sum1=0.0;
  double sum2=0.0;
 
  if(mech == NULL) {
    printf("ERROR: NULL zerork::mechanism pointer in\n");
    printf("       int ConstVolume_CT(...)\n");
    return -1;
  }
  const int num_species = mech->getNumSpecies();

  // Note that d_conc_dt[] first stores the specific heat array
  mech->getCv_R_IntEnergy_RT(temp,&d_conc_dt[0],&workspace[0]);

  // Compute the sum of (species concentration * non-dimensional Cv)
  for(j=0; j<num_species; ) {
    sum1+=conc[j]*d_conc_dt[j];
    ++j;
  }

  mech->getNetReactionRates(temp,&conc[0],&d_conc_dt[0]);
  // Compute the sum of (non-dimensional internal energy *
  //                     species concentration time derivative)
  for(j=0; j<num_species; ) {
    sum2+=workspace[j]*d_conc_dt[j];
    ++j;
  }
  (*d_temp_dt) = -temp*sum2/sum1;
  return 0;
}

// ConstVolume_CMT(...)
//
// Computes the time derivatives of the composition and thermodynamic state
// for a constant volume, homogeneous reactor based on the CMT state
// specification (C = concentration, M = mixture concentration, 
// T = temperature).  Note that the sum of the concentrations is treated 
// independently for the purposes of calculating the pressure depenedent rates.
//
// Inputs (unmodified):
//   zerork::mechanism *mech - mechanism object
//   const double conc[] - species concentration [kmol/m^3], 
//                         length=mech->numSpecies()
//   const double mix_conc - mixture concentration [kmol/m^3]
//   const double temp   - temperature [K]
//  
// Outputs:
//   double workspace[] - workspace array 
//                        min length=3*mech->numSpecies()+mech->numSteps()
//                        Upon return,
//                        workspace[0:numSteps()-1] 
//                            = step rate of progress
//                        workspace[numSteps():numSteps()+numSpecies()-1]
//                            = creation rate of each species
//                        workspace[next numSpecies()]
//                            = destruction rate of each species
//                        workspace[next numSpecies()]
//                            = non-dimensional internal energy of each species
//   double d_conc_dt[] - time derivative of the concentrations [kmol/m^3/s]
//                        length=mech->numSteps()
//   double d_mix_conc_dt - time derivative of the mixture concentration
//                          [kmol/m^3/s]
//   double *d_temp_dt  - pointer to time derivative of temperature [K/s]
//
// Return:
//   int errorFlag      - returns 0 if no error, or -1 if there is an error
int ConstVolume_CMT(zerork::mechanism *mech,
                    const double conc[],
                    const double mix_conc,
                    const double temp,
                    double workspace[],
                    double d_conc_dt[],
                    double *d_mix_conc_dt,
                    double *d_temp_dt)
{
  int j;
  double sum1=0.0;
  double sum2=0.0;
  double *internal_energy,*creation_rate,*destruction_rate, *step_rate;
 
  if(mech == NULL) {
    printf("ERROR: NULL zerork::mechanism pointer in\n");
    printf("       int ConstVolume_CMT(...)\n");
    return -1;
  }
  const int num_species = mech->getNumSpecies();
  const int num_steps = mech->getNumSteps();
  step_rate        = &workspace[0];
  creation_rate    = &workspace[num_steps];
  destruction_rate = &workspace[num_steps+num_species];
  internal_energy  = &workspace[num_steps+2*num_species];

  // Note that d_conc_dt[] first stores the specific heat array
  mech->getCv_R_IntEnergy_RT(temp,&d_conc_dt[0], internal_energy);

  // Compute the sum of (species concentration * non-dimensional Cv)
  for(j=0; j<num_species; ) {
    sum1+=conc[j]*d_conc_dt[j];
    ++j;
  }

  mech->getReactionRatesFromTCM(temp,
                                &conc[0],
                                mix_conc,
                                &d_conc_dt[0],
                                creation_rate,
                                destruction_rate,
                                step_rate);
  // Compute the sum of (non-dimensional internal energy *
  //                     species concentration time derivative)
  for(j=0; j<num_species; ) {
    sum2+=internal_energy[j]*d_conc_dt[j];
    ++j;
  }
  (*d_temp_dt) = -temp*sum2/sum1;

  (*d_mix_conc_dt) = 0.0;
  for(j=0; j<num_species; ) {
    (*d_mix_conc_dt)+=d_conc_dt[j];
    ++j;
  }
  return 0;
}

} // namespace zerork

