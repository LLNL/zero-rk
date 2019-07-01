
#include <stdlib.h>
#include <math.h>

#include <limits>

#include "derivative/const_vol.h"
#include "jacobian_cmt.h"

namespace zerork {

int ConstVolumeSparseJacobian_CMT(mechanism *mechp,
                                  ExplodedROPDerivative *exploded_deriv,
                                  const double conc_orig[],
                                  const double mix_orig,
                                  const double temp_orig,
                                  const double min_conc,
                                  const double state_deriv[],
                                  const int num_non_zeros,
                                  const int column_sum[],
                                  const int row_id[], 
                                  double workspace[],
                                  double jacobian[])
{
  const int num_species   = mechp->getNumSpecies();
  const int num_steps     = mechp->getNumSteps();
  const int min_workspace_size = 9*num_species+num_steps;
  const double sqrt_epsilon = sqrt(numeric_limits<double>::epsilon());

  if(workspace == NULL) {
    return min_workspace_size;
  }
  int num_terms,sparse_id,id; 
  double dM_dt_pos,dM_dt_dtemp,dM_dt_dmix,dM_dt_orig;
  double dT_dt_pos,dT_dt_dtemp,dT_dt_dmix,dT_dt_orig;
  double temp_perturb,dtemp,mix_perturb,dmix,inv_delta;
  double cv_sum,jacobian_sum,mult1,mult2;

  double *conc_pos         = &workspace[0];
  double *inv_conc         = &workspace[num_species];
  double *net_rate_pos     = &workspace[2*num_species];
  double *net_rate_dtemp   = &workspace[3*num_species];
  double *net_rate_dmix    = &workspace[4*num_species];
  double *net_rate_orig    = &workspace[5*num_species];
  // ode_workspace needs a minimum length of 3*numSpecies()+numSteps()
  double *ode_workspace    = &workspace[6*num_species];
  double *specific_heat_cv = &workspace[6*num_species];
  double *internal_energy  = &workspace[7*num_species];

  // compute the strictly positive species concentration
  for(int j=0; j<num_species; ++j) {
    conc_pos[j] = conc_orig[j];
    if(conc_pos[j] < min_conc) {
      conc_pos[j] = min_conc;
    }
    inv_conc[j]   = 1.0/conc_pos[j];
  }
  // compute the species Jacobian from the positive species
  ConstVolume_CMT(mechp,
                  conc_pos,
                  mix_orig,
                  temp_orig,
                  ode_workspace,
                  net_rate_pos,
                  &dM_dt_pos,
                  &dT_dt_pos);

  // step rate of progress is returned in the first numSteps elements of
  // the ode_workspace array
  exploded_deriv->GetSpeciesJacobian(&inv_conc[0],
                                     &ode_workspace[0],
                                     &jacobian[0]);
  // at this point jacobian[] stores the species Jacobian d(wdot[k])/dC[j]
  // ignoring the contribution of perturbations in enhanced third bodies.
 
  // recompute or copy the ODE derivatives (if provided)
  if(state_deriv==NULL) {
    ConstVolume_CMT(mechp,
                    conc_orig,
                    mix_orig,
                    temp_orig,
                    ode_workspace,
                    net_rate_orig,
                    &dM_dt_orig,
                    &dT_dt_orig);
  } else {
    for(int j=0; j<num_species; ++j) {
      net_rate_orig[j] = state_deriv[j];
    }
    dM_dt_orig = state_deriv[num_species];
    dT_dt_orig = state_deriv[num_species+1];
  }
  //-------------------------------------------------------------------------
  // compute the perturbation in concentration
  dmix = mix_orig*sqrt_epsilon;
  mix_perturb = mix_orig + dmix;
  printf("C_mix (orig): %.18g, C_mix (perturb): %.18g\n",
         mix_orig,mix_perturb);
  dmix = mix_perturb - mix_orig;
  // compute the ODE derivative from the original species concentration and
  // temperature using a perturbed mixture concentration 
  ConstVolume_CMT(mechp,
                  conc_orig,
                  mix_perturb,
                  temp_orig,
                  ode_workspace,
                  net_rate_dmix,
                  &dM_dt_dmix,
                  &dT_dt_dmix);
  // update the d/dmix column of the jacobian
  inv_delta = 1.0/dmix;
  sparse_id = column_sum[num_species];
  num_terms = column_sum[num_species+1]-column_sum[num_species];
  for(int j=0; j<num_terms; ++j) {
    id = row_id[sparse_id];
    if(id < num_species) {
      jacobian[sparse_id] = inv_delta*(net_rate_dmix[id]-net_rate_orig[id]);
    } else if(id == num_species) {
      jacobian[sparse_id] = inv_delta*(dM_dt_dmix-dM_dt_orig);
    } else if(id == num_species+1) {
      jacobian[sparse_id] = inv_delta*(dT_dt_dmix-dT_dt_orig);
    } else {
      printf("ERROR: In ConstVolumeSparseJacobian_CMT(...),\n");
      printf("       can not access row %d for sparse index %d\n",
             id,sparse_id);
      exit(-1);
    }
    ++sparse_id;
  }    
  //--------------------------------------------------------------------------
  // compute the perturbation in temperature
  dtemp = temp_orig*sqrt_epsilon;
  temp_perturb = temp_orig + dtemp;
  printf("T [K] (orig): %.18g, T [K] (perturb): %.18g\n",
         temp_orig,temp_perturb);
  dtemp = temp_perturb - temp_orig;
  // compute the ODE derivative from the original species concentration and
  // temperature using a perturbed mixture concentration 
  ConstVolume_CMT(mechp,
                  conc_orig,
                  mix_orig,
                  temp_perturb,
                  ode_workspace,
                  net_rate_dtemp,
                  &dM_dt_dtemp,
                  &dT_dt_dtemp);
  // update the d/dtemp column of the jacobian
  inv_delta = 1.0/dtemp;
  sparse_id = column_sum[num_species+1];
  num_terms = column_sum[num_species+2]-column_sum[num_species+1];
  for(int j=0; j<num_terms; ++j) {
    id = row_id[sparse_id];
    if(id < num_species) {
      jacobian[sparse_id] = inv_delta*(net_rate_dtemp[id]-net_rate_orig[id]);
    } else if(id == num_species) {
      jacobian[sparse_id] = inv_delta*(dM_dt_dtemp-dM_dt_orig);
    } else if(id == num_species+1) {
      jacobian[sparse_id] = inv_delta*(dT_dt_dtemp-dT_dt_orig);
    } else {
      printf("ERROR: In ConstVolumeSparseJacobian_CMT(...),\n");
      printf("       can not access row %d for sparse index %d\n",
             id,sparse_id);
      exit(-1);
    }
    ++sparse_id;
  }
  // This is the point to modify the species jacobian for the 3rd body
  // enhancements
  // :
  // :
  // :
    
  //--------------------------------------------------------------------------
  // compute the d/dconc[j] portion of the row for the mixture concentration
  for(int j=0; j<num_species; ++j) { // column j
    sparse_id = column_sum[j];
    num_terms = column_sum[j+1]-column_sum[j];
    jacobian_sum = 0.0;
    for(int k=0; k<num_terms; ++k) {
      id = row_id[sparse_id];
      if(id < num_species) {
        jacobian_sum += jacobian[sparse_id];
      } else if(id == num_species) {
        jacobian[sparse_id] = jacobian_sum;
      }
      ++sparse_id;
    } // k^th sparse term
  }
  // compute the d/dconc[j] portion of the row for temperature
  
  mechp->getCv_R_IntEnergy_RT(temp_orig,
                              specific_heat_cv,
                              internal_energy);
  cv_sum=0.0;
  for(int j=0; j<num_species; ++j) {
    cv_sum += conc_pos[j]*specific_heat_cv[j];
  }
  mult1 = -temp_orig/cv_sum;
  mult2 = -dT_dt_pos/cv_sum;

  for(int j=0; j<num_species; ++j) { // column j
    sparse_id = column_sum[j];
    num_terms = column_sum[j+1]-column_sum[j];
    jacobian_sum = 0.0;
    for(int k=0; k<num_terms; ++k) {
      id = row_id[sparse_id];
      if(id < num_species) {
        jacobian_sum += jacobian[sparse_id]*internal_energy[id];
      } else if(id == num_species+1) {
        jacobian[sparse_id] = mult1*jacobian_sum + mult2*specific_heat_cv[j];
      }
      ++sparse_id;
    } // k^th sparse term
  }

  return 0;
}

} // namespace zerork
