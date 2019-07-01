#ifndef DERIVATIVE_CONST_VOL_H
#define DERIVATIVE_CONST_VOL_H

#include "zerork/mechanism.h"

namespace zerork {

// constVolume_ct(...)
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
//   double workspace[] - workspace array length=mech->numSpecies()
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
                   double *d_temp_dt);

int ConstVolume_CMT(zerork::mechanism *mech,
                    const double conc[],
                    const double mix_conc,
                    const double temp,
                    double workspace[],
                    double d_conc_dt[],
                    double *d_mix_conc_dt,
                    double *d_temp_dt);
} // namespace zerork

#endif
