#include "kinsol_functions.h"
#include "flame_params.h"

#include <cassert>
#include <algorithm> //min/max

extern "C" void dgbtrf_(int* dim1, int* dim2, int* nu, int* nl, double* a, int* lda, int* ipiv, int* info);
extern "C" void dgbtrs_(char *TRANS, int *N, int *NRHS, int* nu, int* nl, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);


int ConstPressureFlame(N_Vector y,
		       N_Vector ydot, // ydot is the residual
		       void *user_data)
{
  FlameParams *params = (FlameParams *)user_data;
  const int num_local_points = params->num_local_points_;
  const int num_states  = params->reactor_->GetNumStates();
  int Nlocal = num_local_points*num_states;
  int flag = 0;

  // Parallel communications
  flag = ConstPressureFlameComm(Nlocal, y, user_data);
  if(flag != 0) return flag;
  // RHS calculations
  flag = ConstPressureFlameLocal(Nlocal, y, ydot, user_data);

  return flag;
}

int ConstPressureFlameComm(int nlocal,
			   N_Vector y,
			   void *user_data)
{
  FlameParams *params = (FlameParams *)user_data;
  double *y_ptr    = N_VGetArrayPointer(y);
  const int num_local_points = params->num_local_points_;
  const int num_states  = params->reactor_->GetNumStates();
  const int num_species = params->num_species_;
  int my_pe = params->my_pe_;
  int npes  = params->npes_;

#ifdef ZERORK_MPI
  MPI_Comm comm = params->comm_;
  MPI_Status status;
#endif
  const int nover = params->nover_;
  long int dsize = num_states*nover;
  long int dsize_rel = nover;

  // Compute relative volume
  const double RuTref_p = params->reactor_->GetGasConstant()*
    params->reference_temperature_/params->pressure_;

  for(int j=0; j<num_local_points; ++j) {
    int temp_id = j*num_states+num_species + 1;

    double mass_fraction_sum = 0.0;
    for(int k=0; k<num_species; ++k)
      mass_fraction_sum += params->inv_molecular_mass_[k]*y_ptr[j*num_states+k];

    params->rel_vol_[j] = RuTref_p*y_ptr[temp_id]*mass_fraction_sum;
  }

  // Copy y_ptr data into larger arrays
  for (int j=0; j<num_states*num_local_points; ++j)
    params->y_ext_[num_states*nover + j] = y_ptr[j];
  for (int j=0; j<num_local_points; ++j)
    params->rel_vol_ext_[nover+j] = params->rel_vol_[j];

#ifdef ZERORK_MPI
  // MPI sendrecv
  int nodeDest = my_pe-1;
  if (nodeDest < 0) nodeDest = npes-1;
  int nodeFrom = my_pe+1;
  if (nodeFrom > npes-1) nodeFrom = 0;
  MPI_Sendrecv(&params->y_ext_[nover*num_states], dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
	       &params->y_ext_[num_states*(num_local_points+nover)], dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
  MPI_Sendrecv(&params->rel_vol_ext_[nover], dsize_rel, PVEC_REAL_MPI_TYPE, nodeDest, 0,
  	       &params->rel_vol_ext_[num_local_points+nover], dsize_rel, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);

  nodeDest = my_pe+1;
  if (nodeDest > npes-1) nodeDest = 0;
  nodeFrom = my_pe-1;
  if (nodeFrom < 0) nodeFrom = npes-1;
  MPI_Sendrecv(&params->y_ext_[num_states*num_local_points], dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
	       &params->y_ext_[0], dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
  MPI_Sendrecv(&params->rel_vol_ext_[num_local_points], dsize_rel, PVEC_REAL_MPI_TYPE, nodeDest, 0,
  	       &params->rel_vol_ext_[0], dsize_rel, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
#endif

  return 0;
}

int ConstPressureFlameLocal(int nlocal,
			    N_Vector y,
			    N_Vector ydot,
			    void *user_data)
{
  FlameParams *params = (FlameParams *)user_data;
  double *y_ptr    = N_VGetArrayPointer(y);   // caution: assumes realtype == double
  double *ydot_ptr = N_VGetArrayPointer(ydot); // caution: assumes realtype == double

  const int num_local_points = params->num_local_points_;
  const int num_states  = params->reactor_->GetNumStates();
  const int num_species = params->num_species_;
  const int num_local_states = num_local_points*num_states;
  const int convective_scheme_type = params->convective_scheme_type_;

#ifdef ZERORK_MPI
  MPI_Comm comm = params->comm_;
#endif
  int my_pe = params->my_pe_;
  int npes  = params->npes_;

  const double ref_temperature = params->reference_temperature_;

  // Create arrays for RHS, Conv is for convective term
  std::vector<double> rhs, rhsConv;
  rhs.assign(num_local_points*num_states,0.0);
  rhsConv.assign(num_local_points*num_states,0.0);

  double cp_flux_sum, mass_fraction_sum;
  double relative_volume_j;
  int transport_error;

  double local_sum;
  double local_max;

  double velocity, thermal_diffusivity;

  int temp_out_of_bounds = 0;
  int global_temp_out_of_bounds = 0;
  for(int j=0; j<num_local_points; ++j) {
    double temperature = y_ptr[j*num_states + num_states-1]*params->reference_temperature_;
    if(temperature < 100.0 || temperature > 10000) {
      temp_out_of_bounds += 1;
    }
  }
#ifdef ZERORK_MPI
  MPI_Allreduce(&temp_out_of_bounds,&global_temp_out_of_bounds,1,MPI_INT,MPI_MAX,comm);
#else
  global_temp_out_of_bounds = temp_out_of_bounds;
#endif
  if(global_temp_out_of_bounds > 0) return 1;// recoverable error

  // Set the residual and rhs to zero
  for(int j=0; j<num_local_states; ++j) {
    ydot_ptr[j] = 0.0;
    rhs[j] = 0.0;
    rhsConv[j] = 0.0;
  }

  // Compute relative volume
  const double RuTref_p = params->reactor_->GetGasConstant()*
    params->reference_temperature_/params->pressure_;

  for(int j=0; j<num_local_points; ++j) {
    int temp_id = j*num_states+num_species + 1;

    mass_fraction_sum = 0.0;
    for(int k=0; k<num_species; ++k)
      mass_fraction_sum += params->inv_molecular_mass_[k]*y_ptr[j*num_states+k];

    params->rel_vol_[j] = RuTref_p*y_ptr[temp_id]*mass_fraction_sum;
  }

  //--------------------------------------------------------------------------
  // Copy data in extended work arrays with ghost cells
  const int nover = params->nover_;

  // Larger arrays with ghost cells
  std::vector<double> dz, dzm, inv_dz, inv_dzm;
  dz.assign( num_local_points+(2*nover), 0.0);
  dzm.assign( num_local_points+(2*nover), 0.0);
  inv_dz.assign( num_local_points+(2*nover), 0.0);
  inv_dzm.assign( num_local_points+(2*nover), 0.0);

  // Copy y_ptr data into larger arrays
  for (int j=0; j<num_states*num_local_points; ++j)
    params->y_ext_[num_states*nover + j] = y_ptr[j];

  for (int j=0; j<num_local_points; ++j)
    params->rel_vol_ext_[nover+j] = params->rel_vol_[j];

  for (int j=0; j<num_local_points+2*nover; ++j) {
    dz[j] = params->dz_local_[j];
    dzm[j] = params->dzm_local_[j];
    inv_dz[j] = params->inv_dz_local_[j];
    inv_dzm[j] = params->inv_dzm_local_[j];
  }

  //--------------------------------------------------------------------------
  // Set boundary conditions
  // First proc: inlet conditions in ghost cells
  if (my_pe == 0) {
    for(int j=0; j<nover; ++j) {
      for(int k=0; k<num_species; ++k) {
	params->y_ext_[j*num_states + k] = params->inlet_mass_fractions_[k];
      }
      params->y_ext_[(j+1)*num_states-1] = params->inlet_temperature_;
      params->y_ext_[(j+1)*num_states-2] = y_ptr[num_states-2];
      params->rel_vol_ext_[j] = params->inlet_relative_volume_;
    }
  }

  // Last proc: zero gradient. Copy last domain cell in ghost cells
  if (my_pe == npes-1) {
    for(int j=num_local_points+nover; j<num_local_points+2*nover; ++j) {
      for(int k=0; k<num_species; ++k) {
	params->y_ext_[j*num_states + k] = y_ptr[(num_local_points-1)*num_states+k];
      }
      params->y_ext_[(j+1)*num_states-1] = y_ptr[num_local_points*num_states-1];
      params->y_ext_[(j+1)*num_states-2] = y_ptr[num_local_points*num_states-2];
      params->rel_vol_ext_[j] = params->rel_vol_[num_local_points-1];
    }
  }

  //--------------------------------------------------------------------------
  // Compute the constant pressure reactor source term
  for(int j=0; j<num_local_points; ++j) {

    params->reactor_->GetTimeDerivativeSteady(&y_ptr[j*num_states],
                                              &params->step_limiter_[0],
                                              &rhs[j*num_states]);
  }

  //--------------------------------------------------------------------------
  // Compute the interior heat capacity, conductivity, and species mass fluxes.
  for(int j=0; j<num_local_points+1; ++j) {
    int jext = j + nover;

    // compute the upstream mid point state for the transport calculations
    for(int k=0; k<num_species; ++k) {

      // mid point mass fractions
      params->transport_input_.mass_fraction_[k] =
	0.5*(params->y_ext_[jext*num_states+k] + params->y_ext_[(jext-1)*num_states+k]);

      // mid point mass fraction gradient
      params->transport_input_.grad_mass_fraction_[k] = inv_dz[jext]*
	(params->y_ext_[jext*num_states+k] - params->y_ext_[(jext-1)*num_states+k]);
    }

    // mid point temperature
    params->transport_input_.temperature_ = 0.5*ref_temperature*
      (params->y_ext_[(jext+1)*num_states-1] + params->y_ext_[jext*num_states-1]);

    // mid point temperature gradient
    params->transport_input_.grad_temperature_[0] = inv_dz[jext]*ref_temperature*
      (params->y_ext_[(jext+1)*num_states-1] - params->y_ext_[jext*num_states-1]);

    // mixture specific heat at mid point. Species cp will be overwritten
    // for frozen thermo only
    params->mixture_specific_heat_mid_[j] =
      params->reactor_->GetMixtureSpecificHeat_Cp(
		       		  params->transport_input_.temperature_,
				  &params->transport_input_.mass_fraction_[0],
				  &params->species_specific_heats_[0]);

    // Reset species cp
    for(int k=0; k<num_species; k++) {
      params->species_specific_heats_[num_species*j+k] = 0.0;
    }

    //mixture molecular mass at mid point
    // for frozen thermo only
    double mass_fraction_weight_sum = 0.0;
    for(int k=0; k<num_species; ++k) {
      mass_fraction_weight_sum += params->inv_molecular_mass_[k]*params->transport_input_.mass_fraction_[k];
    }
    params->molecular_mass_mix_mid_[j] = 1.0/mass_fraction_weight_sum;

    // compute the conductivity at the upstream mid point (j-1/2)
    transport_error = params->transport_->GetMixtureConductivity(
	      			     params->transport_input_,
				     &params->thermal_conductivity_[j]);
    if(transport_error != transport::NO_ERROR) {
      return transport_error;
    }

    // specific heat at grid point j
    if (j != num_local_points) { //not used in derivatives
      params->mixture_specific_heat_[j] =
	params->reactor_->GetMixtureSpecificHeat_Cp(
				  ref_temperature*params->y_ext_[(jext+1)*num_states-1],
				  &params->y_ext_[jext*num_states],
				  &params->species_specific_heats_[num_species*j]);
    }
    /*
    // compute the species mass flux at the upstream mid point
    transport_error = params->transport_->GetSpeciesMassFlux(
				   params->transport_input_,
				   num_species,
				   &params->thermal_conductivity_[j],
				   &params->mixture_specific_heat_mid_[j],
				   &params->species_mass_flux_[j*num_species],
				   &params->species_lewis_numbers_[j*num_species]);
    */
    // Frozen thermo is the default
    // TODO: switch between regular GetSpeciesMassFlux and FrozenThermo from the input file
    /**/
    transport_error = params->transport_->GetSpeciesMassFluxFrozenThermo(
					 params->transport_input_,
					 num_species,
					 &params->thermal_conductivity_[j],
					 &params->mixture_specific_heat_mid_[j],
					 &params->species_mass_flux_[j*num_species],
					 &params->species_lewis_numbers_[j*num_species]);
    /**/

    if(transport_error != transport::NO_ERROR) {
      return transport_error;
    }

  } // for j<num_local_points+1

  //--------------------------------------------------------------------------
  // Compute convective and diffusive terms for species and temperature
  for(int j=0; j<num_local_points; ++j) {
    int jext = j + nover;
    int jglobal = j + my_pe*num_local_points;

    double b=0,c=0,d=0,e=0; //coefficients of j+1, j, j-1, j-2 terms
    assert(convective_scheme_type >= 0 && convective_scheme_type < 3);
    if(convective_scheme_type == 0) {
      // First order upwind
      b = 0;
      c =  inv_dz[jext];
      d = -inv_dz[jext];
      e = 0;
    } else if(convective_scheme_type == 1) {
      // Second order upwind
      b = 0;
      c = inv_dz[jext] + 1.0/(dz[jext]+dz[jext-1]);
      d = -(dz[jext]+dz[jext-1])/(dz[jext]*dz[jext-1]);
      e = dz[jext]/dz[jext-1]/(dz[jext]+dz[jext-1]);
    } else if(convective_scheme_type == 2) {
      // Centered
      b = dz[jext]/dz[jext+1]/(dz[jext]+dz[jext+1]);
      c = (dz[jext+1]-dz[jext])/dz[jext+1]/dz[jext];
      d = -dz[jext+1]/dz[jext]/(dz[jext]+dz[jext+1]);
      e = 0;
    }

    relative_volume_j   = params->rel_vol_ext_[jext];

    // compute the species mass fraction advection and diffusion
    for(int k=0; k<num_species; ++k) {
      rhsConv[j*num_states+k] -= relative_volume_j*
	(b*params->y_ext_[(jext+1)*num_states+k] +
	 c*params->y_ext_[jext*num_states+k] +
	 d*params->y_ext_[(jext-1)*num_states+k]);

      rhs[j*num_states+k] -= relative_volume_j*inv_dzm[jext]*
	( params->species_mass_flux_[num_species*(j+1)+k]
	  -params->species_mass_flux_[num_species*j+k]);
    }

    // compute the species specific heat diffusive flux sum
    cp_flux_sum = 0.0;
    for(int k=0; k<num_species; ++k) {
      cp_flux_sum += params->species_specific_heats_[num_species*j+k]*
	0.5*(params->species_mass_flux_[num_species*j+k]+
	     params->species_mass_flux_[num_species*(j+1)+k]);
    }

    // Compute the temperature advection
    rhsConv[(j+1)*num_states-1] -= relative_volume_j*
      (b*params->y_ext_[(jext+2)*num_states-1] +
       c*params->y_ext_[(jext+1)*num_states-1] +
       d*params->y_ext_[jext*num_states-1]);

    // Flux term
    rhs[(j+1)*num_states-1] -= relative_volume_j*
      cp_flux_sum/params->mixture_specific_heat_[j]*
      (b*params->y_ext_[(jext+2)*num_states-1] +
       c*params->y_ext_[(jext+1)*num_states-1] +
       d*params->y_ext_[jext*num_states-1]);

    // Add the thermal diffusivity contribution to dT[j]/dt
    rhs[(j+1)*num_states-1] +=
      (relative_volume_j*inv_dzm[jext]/params->mixture_specific_heat_[j])*
      (params->thermal_conductivity_[j+1]*inv_dz[jext+1]*
       (params->y_ext_[(jext+2)*num_states-1] - params->y_ext_[(jext+1)*num_states-1])
       -params->thermal_conductivity_[j]*inv_dz[jext]*
       (params->y_ext_[(jext+1)*num_states-1] - params->y_ext_[jext*num_states-1]));

    // Compute the mass flow derivative
    rhs[(j+1)*num_states-2] = 0.0;
    if (jglobal == params->j_fix_) {
      rhs[(j+1)*num_states-2] = (y_ptr[(j+1)*num_states-1] - params->temperature_fix_);
    }
    else if (jglobal > params->j_fix_) {
      rhs[(j+1)*num_states-2] = (params->y_ext_[(jext+1)*num_states-2] -
				 params->y_ext_[jext*num_states-2])*inv_dz[jext];
    }
    else if (jglobal < params->j_fix_) {
      rhs[(j+1)*num_states-2] = -(params->y_ext_[(jext+2)*num_states-2] -
				 params->y_ext_[(jext+1)*num_states-2])*inv_dz[jext+1];
    }

  } // for(int j=1; j<num_local_points; ++j) // loop computing rhs

  // -------------------------------------------------------------------------
  // Compute the final residuals
  for(int j=0; j<num_local_points; ++j) {
    int mflux_id = j*num_states+num_species;
    int temp_id = mflux_id + 1;

    for(int k=0; k<num_species; ++k)
      ydot_ptr[j*num_states+k] = rhsConv[j*num_states+k]*y_ptr[mflux_id] + rhs[j*num_states+k];

    ydot_ptr[mflux_id] = rhs[mflux_id];
    ydot_ptr[temp_id] = rhsConv[temp_id]*y_ptr[mflux_id] + rhs[temp_id];

    // Copy rhsConv in params for use in jacobian
    for(int k=0; k<num_states; ++k)
      params->rhsConv_[j*num_states + k] = rhsConv[j*num_states+k];

    // Switch RHS at j_fix to have non-zero diagonal
    int jglobal = j + my_pe*num_local_points;
    if (jglobal == params->j_fix_) {
      ydot_ptr[mflux_id] = ydot_ptr[temp_id];
      ydot_ptr[temp_id] = rhs[mflux_id];
    }

  } // for j<num_local_points

  // Add time derivative term if pseudo unsteady
  if(params->pseudo_unsteady_) {
    for(int j=0; j<num_local_points; ++j) {
      int mflux_id  = j*num_states+num_species; // relative volume index of pt j
      int temp_id   = mflux_id+1;               // temperature index of pt j

      for(int k=0; k<num_species; ++k) {
        ydot_ptr[j*num_states+k] -= (y_ptr[j*num_states+k] - params->y_old_[j*num_states+k])/
          params->dt_;
      }

      ydot_ptr[temp_id] -= (y_ptr[temp_id] - params->y_old_[temp_id])/params->dt_;
    }
  }

  //------------------------------------------------------------------
  // Parallel communication for finite difference jacobian
  // TO DO: Move to a separate function?
  if(params->integrator_type_ == 2 || params->integrator_type_ == 3) {
#ifdef ZERORK_MPI
    MPI_Status status;
#endif
    long int dsize = num_states*nover;

    // Copy y_ptr and ydot_ptr into larger arrays
    for (int j=0; j<num_states*num_local_points; ++j)
      params->rhs_ext_[num_states*nover + j] = ydot_ptr[j];
#ifdef ZERORK_MPI
    // MPI sendrecv
    int nodeDest = my_pe-1;
    if (nodeDest < 0) nodeDest = npes-1;
    int nodeFrom = my_pe+1;
    if (nodeFrom > npes-1) nodeFrom = 0;
    MPI_Sendrecv(&params->rhs_ext_[nover*num_states], dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
		 &params->rhs_ext_[num_states*(num_local_points+nover)], dsize, PVEC_REAL_MPI_TYPE,
		 nodeFrom, 0, comm, &status);

    nodeDest = my_pe+1;
    if (nodeDest > npes-1) nodeDest = 0;
    nodeFrom = my_pe-1;
    if (nodeFrom < 0) nodeFrom = npes-1;
    MPI_Sendrecv(&params->rhs_ext_[num_states*num_local_points], dsize, PVEC_REAL_MPI_TYPE, nodeDest,
		 0, &params->rhs_ext_[0], dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
#endif
  }

  // -------------------------------------------------------------------------
  // Post-processing for output
  // Compute laminar flame speed from mass flux (only on rank = 0)
  if(my_pe == 0) {
    params->flame_speed_ = y_ptr[num_species]*params->inlet_relative_volume_;
  }
#ifdef ZERORK_MPI
  MPI_Bcast(&params->flame_speed_, 1, MPI_DOUBLE, 0, comm);
#endif

  // Compute flame thick l_F = (T_max - T_min)/|gradT|_max
  double local_temperature;
  local_max = 0.0;
  for(int j=0; j<num_local_points; ++j) {
    int jext = j + nover;
    local_temperature = ref_temperature*params->y_ext_[(jext+1)*num_states-1];
    if(local_temperature > local_max) {
      local_max = local_temperature;
    }
  }
#ifdef ZERORK_MPI
  MPI_Allreduce(&local_max,&params->max_temperature_,1,PVEC_REAL_MPI_TYPE,MPI_MAX,comm);
#else
  params->max_temperature_ = local_max;
#endif
  double gradT;
  local_max = 0.0;
  for(int j=0; j<num_local_points; ++j) {
    int jext = j + nover;
    gradT = fabs(inv_dz[jext]*ref_temperature*
		 (params->y_ext_[(jext+1)*num_states-1] - params->y_ext_[jext*num_states-1]));
    if (gradT > local_max) {
      local_max = gradT;
    }
  }
#ifdef ZERORK_MPI
  MPI_Allreduce(&local_max,&params->flame_thickness_,1,PVEC_REAL_MPI_TYPE,MPI_MAX,comm);
#else
  params->flame_thickness_ = local_max;
#endif
  params->flame_thickness_ = (params->max_temperature_-params->inlet_temperature_)/
                                 params->flame_thickness_;

  // Compute alternate flame thickness definition lF = K/rho/cp/SL (alpha = K/rho/cp)
  // using inlet values
  if(my_pe == 0) {
    params->flame_thickness_alpha_ = params->thermal_conductivity_[0]/
      params->inlet_relative_volume_/params->mixture_specific_heat_[0]/params->flame_speed_;
  }
#ifdef ZERORK_MPI
  MPI_Bcast(&params->flame_thickness_alpha_, 1, MPI_DOUBLE, 0, comm);
#endif

  return 0;
}

//------------------------------------------------------------------
// Banded Block Diagonal preconditioner, factorized with SuperLU
#if defined SUNDIALS2
int ReactorBBDSetup(N_Vector y, // [in] state vector
		    N_Vector yscale, // [in] state scaler
		    N_Vector ydot, // [in] state derivative
		    N_Vector ydotscale, // [in] state derivative scaler
		    void *user_data, // [in/out]
                    N_Vector tmp1, N_Vector tmp2)
{
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorBBDSetup(N_Vector y, // [in] state vector
		    N_Vector yscale, // [in] state scaler
		    N_Vector ydot, // [in] state derivative
		    N_Vector ydotscale, // [in] state derivative scaler
		    void *user_data) // [in/out]
{
#endif
  FlameParams *params    = (FlameParams *)user_data;
  const int num_local_points   = params->num_local_points_;
  const int num_states   = params->reactor_->GetNumStates();
  const int num_nonzeros_loc = params->num_nonzeros_loc_;
  const int num_local_states = num_states*num_local_points;
  const int num_total_points = params->num_points_;
  const int num_total_states = num_states*num_total_points;
  double *y_ptr          = N_VGetArrayPointer(y);
  double *ydot_ptr       = N_VGetArrayPointer(ydot);
  int error_flag = 0;
  double alpha = 1.0e-6;
  double beta = 1.0e-14;
  double delta;

  const int my_pe = params->my_pe_;
  const int npes  = params->npes_;
  const int nover=params->nover_;

  // Create work arrays
  std::vector<double> y_saved,rhs_ext_saved;
  y_saved.assign(num_local_points*num_states,0.0);
  rhs_ext_saved.assign((num_local_points+2*nover)*num_states,0.0);

  int group, width;
  int mkeep = params->num_off_diagonals_;
  width = 2*mkeep + 1;

  std::vector<double> jac_bnd;
  jac_bnd.assign( (num_local_points+2*nover)*num_states*width, 0.0);

  // Compute RHS
  ConstPressureFlame(y, ydot, user_data);

  // Save copy of state vector and rhs
  for (int j=0; j<num_local_states; ++j)
    y_saved[j] = y_ptr[j];

  for (int j=0; j<num_states*(num_local_points+2*nover); ++j)
    rhs_ext_saved[j] = params->rhs_ext_[j];

  // Banded jacobian
  int j, jglobal, i1global, i2global, iloc, jext, iext, jstate;
  for (group = 1; group <= width; group++) {

    // Perturb y
    for (jglobal=group-1; jglobal<num_total_states; jglobal+=width) {
      j = jglobal - my_pe*num_local_points*num_states;
      if (j>=0 && j<num_local_states) {
	delta = alpha*y_saved[j] + beta;
	y_ptr[j] = y_saved[j] + delta;
      }
    }//for j=group-1 +=width

    // Compute RHS
    ConstPressureFlame(y, ydot, user_data);

    // Compute jacobian
    // here j is the COLUMN and i is the ROW
    for (jglobal=group-1; jglobal<num_total_states; jglobal+=width) {
      j = jglobal - my_pe*num_local_states;
      jstate = jglobal % num_states;
      if (j>=0 && j<num_local_states) {
	i1global = std::max(0, jglobal - jstate - num_states);
        i2global = std::min(jglobal + (num_states-1 - jstate) + num_states, num_total_states-1);
	jext = j + nover*num_states;
	for (int i=i1global; i<=i2global; i++) {
	  iloc = i-my_pe*num_local_states;
	  iext = iloc + nover*num_states;
	  jac_bnd[jext*width + i-jglobal+mkeep] =
	    (params->rhs_ext_[iext] - rhs_ext_saved[iext]) / (y_ptr[j]- y_saved[j]);
	} // for i < i2
	y_ptr[j] = y_saved[j];
      } //if j exists on this processor
    } // for j=group-1 +=width

  } // for group < width

  // Restore the state and rhs vectors back to original values
  for (int j=0; j<num_local_states; ++j) {
    jext = j + nover*num_states;
    y_ptr[j] = y_saved[j];
    ydot_ptr[j] = rhs_ext_saved[jext];
  }

#ifdef ZERORK_MPI
  // Perform parallel communication of jacobian
  MPI_Comm comm = params->comm_;
  MPI_Status status;
  long int dsize_jac_bnd = nover*num_states*width;

  // MPI sendrecv
  int nodeDest = my_pe-1;
  if (nodeDest < 0) nodeDest = npes-1;
  int nodeFrom = my_pe+1;
  if (nodeFrom > npes-1) nodeFrom = 0;
  MPI_Sendrecv(&jac_bnd[nover*num_states*width], dsize_jac_bnd, PVEC_REAL_MPI_TYPE, nodeDest, 0, &jac_bnd[num_states*(num_local_points+nover)*width], dsize_jac_bnd, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);

  nodeDest = my_pe+1;
  if (nodeDest > npes-1) nodeDest = 0;
  nodeFrom = my_pe-1;
  if (nodeFrom < 0) nodeFrom = npes-1;
  MPI_Sendrecv(&jac_bnd[num_states*num_local_points*width], dsize_jac_bnd, PVEC_REAL_MPI_TYPE, nodeDest, 0, &jac_bnd[0], dsize_jac_bnd, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
#endif

  // Get pattern "manually" for now
  // TODO: find a cleaner way
  int innz=0;
  // here j is the ROW and i is the COLUMN
  for (j=0; j<num_local_states; j++) {
    jglobal = j + my_pe*num_local_states;
    jext = j + nover*num_states;
    int jstate = jglobal % num_states;
    i1global = std::max(0, jglobal - jstate - num_states);
    i2global = std::min(jglobal + (num_states-1 - jstate) + num_states, num_total_states-1);
    for (int i=i1global; i<=i2global; i++) {
      iloc = i-my_pe*num_local_states;
      iext = iloc + nover*num_states;
      int istate = i % num_states;
      int dense_id = num_states*istate + jstate; //i is column and j is row
      //Diagonal block.
      if (i>= jglobal-jstate && i<=jglobal+num_states-1-jstate) {
	if (params->dense_to_sparse_[dense_id] == 1) {
	  params->reactor_jacobian_dist_[innz] = jac_bnd[iext*width + jglobal-i+mkeep];
	  innz++;
	}
      }
      //Off-diagonal blocks
      if (i<jglobal-jstate || i>jglobal+num_states-1-jstate) {
	if (params->dense_to_sparse_offdiag_[dense_id] == 1) {
	  params->reactor_jacobian_dist_[innz] = jac_bnd[iext*width + jglobal-i+mkeep];
	  innz++;
	}
      }
    } // for i1 to i2
  } // for j < num_local_states

  // Factorize with SuperLU (parallel is default, serial if specified in input)
  if(params->superlu_serial_) {
    if(params->sparse_matrix_->IsFirstFactor()) {
      error_flag =
	params->sparse_matrix_->FactorNewPatternCRS(num_nonzeros_loc,
						    &params->col_id_[0],
						    &params->row_sum_[0],
						    &params->reactor_jacobian_dist_[0]);
    } else {
      error_flag =
	params->sparse_matrix_->FactorSamePattern(&params->reactor_jacobian_dist_[0]);
    } //if first factor
#ifdef ZERORK_MPI
  } else {
    if(params->sparse_matrix_dist_->IsFirstFactor_dist()) {
      error_flag =
	params->sparse_matrix_dist_->FactorNewPatternCRS_dist(num_nonzeros_loc,
							      &params->col_id_[0],
							      &params->row_sum_[0],
							      &params->reactor_jacobian_dist_[0]);
    } else {
      error_flag =
	params->sparse_matrix_dist_->FactorSamePatternCRS_dist(num_nonzeros_loc,
							       &params->col_id_[0],
							       &params->row_sum_[0],
							       &params->reactor_jacobian_dist_[0]);
    } //if first factor
#endif
  } // if superlu serial

  return error_flag;
}

// Banded block diagonal finite difference Jacobian, solved with SuperLU
#if defined SUNDIALS2
int ReactorBBDSolve(N_Vector y, // [in] state vector
		    N_Vector yscale, // [in] state scaler
		    N_Vector ydot, // [in] state derivative
		    N_Vector ydotscale, // [in] state derivative scaler
		    N_Vector vv, // [in/out] rhs/solution vector
		    void *user_data, // [in/out]
                    N_Vector tmp)
{
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorBBDSolve(N_Vector y, // [in] state vector
		    N_Vector yscale, // [in] state scaler
		    N_Vector ydot, // [in] state derivative
		    N_Vector ydotscale, // [in] state derivative scaler
		    N_Vector vv, // [in/out] rhs/solution vector
		    void *user_data) // [in/out]
{
#endif
  FlameParams *params = (FlameParams *)user_data;
  double *solution = N_VGetArrayPointer(vv);
  int error_flag = 0;

  if(params->superlu_serial_) {
    error_flag = params->sparse_matrix_->Solve(&solution[0],&solution[0]);
#ifdef ZERORK_MPI
  } else {
    error_flag = params->sparse_matrix_dist_->Solve_dist(&solution[0],&solution[0]);
#endif
  }

  return error_flag;
}

//------------------------------------------------------------------
// Approximate factorization preconditioner
// Chemistry Jacobian is (sparse) block diagonal. Each n_sp x n_sp block
// factorized separately with SuperLU
// Transport is tridiagonal over whole domain, factorized with LAPACK
#if defined SUNDIALS2
int ReactorAFSetup(N_Vector y, // [in] state vector
		   N_Vector yscale, // [in] state scaler
		   N_Vector ydot, // [in] state derivative
		   N_Vector ydotscale, // [in] state derivative scaler
		   void *user_data, // [in/out]
                   N_Vector tmp1, N_Vector tmp2)
{
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorAFSetup(N_Vector y, // [in] state vector
		   N_Vector yscale, // [in] state scaler
		   N_Vector ydot, // [in] state derivative
		   N_Vector ydotscale, // [in] state derivative scaler
		   void *user_data) // [in/out]
{
#endif
  FlameParams *params    = (FlameParams *)user_data;
  const int num_local_points   = params->num_local_points_;
  const int num_states   = params->reactor_->GetNumStates();
  const int num_species = num_states - 2;
  const int num_total_points = params->num_points_;
  const int num_nonzeros_zerod = params->reactor_->GetJacobianSize();
  const int num_states_local = params->num_states_local_;
  double *y_ptr          = N_VGetArrayPointer(y);
  int error_flag = 0;
  bool Tfix = false;
  double constant = 1.0e6;//1.0e6
  int my_pe  = params->my_pe_;
  const int nover = params->nover_;

  // Initialize transport Jacobian
  for(int j=0; j<num_local_points*5*num_states; j++)
    params->banded_jacobian_[j] = 0.0;

  // Get grid spacing
  std::vector<double> dz, dzm, inv_dz, inv_dzm;
  dz.assign( num_local_points+(2*nover), 0.0);
  dzm.assign( num_local_points+(2*nover), 0.0);
  inv_dz.assign( num_local_points+(2*nover), 0.0);
  inv_dzm.assign( num_local_points+(2*nover), 0.0);
  for (int j=0; j<num_local_points+2*nover; ++j) {
    dz[j] = params->dz_local_[j];
    dzm[j] = params->dzm_local_[j];
    inv_dz[j] = params->inv_dz_local_[j];
    inv_dzm[j] = params->inv_dzm_local_[j];
  }

  // Evaluate analytic transport J
  const int convective_scheme_type = params->convective_scheme_type_;
  double relative_volume_j, relative_volume_jp, relative_volume_jm;
  for(int j=0; j<num_local_points; j++) {
    int jext = j + nover;
    int jglobal = j + my_pe*num_local_points;

    double bm=0,c=0,dp=0; //coefficients of  j+1, j, j-1 terms
    assert(convective_scheme_type >= 0 && convective_scheme_type < 3);
    if(convective_scheme_type == 0) {
      // First order upwind
      c =  inv_dz[jext];
    } else if(convective_scheme_type == 1) {
      // Second order upwind
      if (y_ptr[j*num_states+num_species] >= 0) {
	// Use points upstream
	c = inv_dz[jext] + 1.0/(dz[jext]+dz[jext-1]);
      } else {
	// Use points downstream
	c = -inv_dz[jext+1] - 1.0/(dz[jext+1]+dz[jext+2]);
      }
    } else if(convective_scheme_type == 2) {
      // Centered
      c = (dz[jext+1]-dz[jext])/dz[jext+1]/dz[jext];
      bm = dz[jext-1]/dz[jext]/(dz[jext-1]+dz[jext]);
      dp = -dz[jext+2]/dz[jext+1]/(dz[jext+1]+dz[jext+2]);
    }

    relative_volume_j  = params->rel_vol_ext_[jext];
    relative_volume_jp = params->rel_vol_ext_[jext+1];
    relative_volume_jm = params->rel_vol_ext_[jext-1];

    // Species
    for(int k=0; k<num_species; k++) {
      // Diagonal
      params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] =
	(-params->thermal_conductivity_[j+1]*inv_dz[jext+1]/
	 params->mixture_specific_heat_mid_[j+1]/
	 params->species_lewis_numbers_[k] -
	 params->thermal_conductivity_[j]*inv_dz[jext]/
	 params->mixture_specific_heat_mid_[j]/
	 params->species_lewis_numbers_[k])*
	relative_volume_j*inv_dzm[jext] -
	c*y_ptr[j*num_states+num_species]*relative_volume_j;

      // drhs_j-1/dY_j
      if(jglobal > 0) {
	params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 - 1] =
	  (params->thermal_conductivity_[j]*inv_dz[jext]/
	   params->mixture_specific_heat_mid_[j]/
	   params->species_lewis_numbers_[k])*
	  relative_volume_jm*inv_dzm[jext-1]-
	  bm*params->y_ext_[(jext-1)*num_states+num_species]*relative_volume_jm;
      }

      // drhs_j+1/dY_j
      if(jglobal < num_total_points-1) {
	params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 1] =
	  (params->thermal_conductivity_[j+1]*inv_dz[jext+1]/
	   params->mixture_specific_heat_mid_[j+1]/
	   params->species_lewis_numbers_[k])*
	  relative_volume_jp*inv_dzm[jext+1] -
	  dp*params->y_ext_[(jext+1)*num_states+num_species]*relative_volume_jp;
      }

    } // for k<num_species

    // Temperature
    // Diagonal
    params->banded_jacobian_[(num_species+1)*(num_local_points*5) + j*5 + 1 + 2 + 0] =
      (-params->thermal_conductivity_[j+1]*inv_dz[jext+1]/
       params->mixture_specific_heat_mid_[j+1] -
       params->thermal_conductivity_[j]*inv_dz[jext]/
       params->mixture_specific_heat_mid_[j])*
      relative_volume_j*inv_dzm[jext] -
      c*y_ptr[j*num_states+num_species]*relative_volume_j;

    // Fixed temperature point
    if (jglobal == params->j_fix_)
      params->banded_jacobian_[(num_species+1)*(num_local_points*5) + j*5 + 1 + 2 + 0] = 1.0;

    // drhs_j-1/dY_j
    if(jglobal > 0) {
      params->banded_jacobian_[(num_species+1)*(num_local_points*5) + j*5 + 1 + 2 - 1] =
	(params->thermal_conductivity_[j]*inv_dz[jext]/
	 params->mixture_specific_heat_mid_[j])*
	relative_volume_jm*inv_dzm[jext-1] -
	bm*params->y_ext_[(jext-1)*num_states+num_species]*relative_volume_jm;
    }
    if (jglobal == params->j_fix_ + 1)
      params->banded_jacobian_[(num_species+1)*(num_local_points*5) + j*5 + 1 + 2 - 1] = 0.0;

    // drhs_j+1/dY_j
    if(jglobal < num_total_points-1) {
      params->banded_jacobian_[(num_species+1)*(num_local_points*5) + j*5 + 1 + 2 + 1] =
	(params->thermal_conductivity_[j+1]*inv_dz[jext+1]/
	 params->mixture_specific_heat_mid_[j+1])*
	relative_volume_jp*inv_dzm[jext+1] -
	dp*params->y_ext_[(jext+1)*num_states+num_species]*relative_volume_jp;
    }
    if (jglobal == params->j_fix_ - 1)
      params->banded_jacobian_[(num_species+1)*(num_local_points*5) + j*5 + 1 + 2 + 1] = 0.0;


    // Mass flux
    jglobal = j + my_pe*num_local_points;
    if (jglobal == params->j_fix_) {
      params->banded_jacobian_[num_species*(num_local_points*5) + j*5 + 1 + 2 + 0] =
        params->rhsConv_[j*num_states+num_species+1]; //diag
      params->banded_jacobian_[num_species*(num_local_points*5) + j*5 + 1 + 2 + 1] = -inv_dz[jext+1];
      params->banded_jacobian_[num_species*(num_local_points*5) + j*5 + 1 + 2 - 1] = -inv_dz[jext];
    }
    if (jglobal > params->j_fix_) {
      params->banded_jacobian_[num_species*(num_local_points*5) + j*5 + 1 + 2 + 0] = inv_dz[jext];
      params->banded_jacobian_[num_species*(num_local_points*5) + j*5 + 1 + 2 + 1] = -inv_dz[jext+1];
      params->banded_jacobian_[num_species*(num_local_points*5) + j*5 + 1 + 2 - 1] = 0.0;
    }

    if (jglobal < params->j_fix_) {
      params->banded_jacobian_[num_species*(num_local_points*5) + j*5 + 1 + 2 + 0] = inv_dz[jext+1];
      params->banded_jacobian_[num_species*(num_local_points*5) + j*5 + 1 + 2 + 1] = 0.0;
      params->banded_jacobian_[num_species*(num_local_points*5) + j*5 + 1 + 2 - 1] = -inv_dz[jext];
    }

    if (jglobal == params->j_fix_ + 1)
      params->banded_jacobian_[num_species*(num_local_points*5) + j*5 + 1 + 2 - 1] = 0.0;
    if (jglobal == params->j_fix_ - 1)
      params->banded_jacobian_[num_species*(num_local_points*5) + j*5 + 1 + 2 + 1] = 0.0;
  } // for j<num_local_points

  // Local chemistry Jacobian (and mass flux)
  if(params->store_jacobian_) {
    params->saved_jacobian_chem_.assign(num_nonzeros_zerod*num_local_points, 0.0);
    for(int j=0; j<num_local_points; ++j) {
      int jglobal = j + my_pe*num_local_points;
      if(jglobal == params->j_fix_) {
	Tfix = true;
      } else {
	Tfix = false;
      }
      // Get Jacobian
      params->reactor_->GetJacobianSteady(&y_ptr[j*num_states],
					  &params->rhsConv_[j*num_states],
					  Tfix,
                                          &params->step_limiter_[0],
					  &params->saved_jacobian_chem_[j*num_nonzeros_zerod]);

      // dmdot_j/dT_j at j=j_fix
      int jext = j + nover;
      if(jglobal == params->j_fix_) {
	for(int k=0; k<num_nonzeros_zerod; ++k) {
	  int col_id = params->column_id_chem_[k];
	  int row_id = params->row_id_chem_[k];
	  if (col_id==num_species+1 && row_id==num_species) {
	    double c = (dz[jext+1]-dz[jext])/dz[jext+1]/dz[jext];

	    params->saved_jacobian_chem_[j*num_nonzeros_zerod+k] =
	      params->saved_jacobian_chem_[j*num_nonzeros_zerod+params->diagonal_id_chem_[num_states-1]] +
	      (-params->thermal_conductivity_[j+1]*inv_dz[jext+1]/
	       params->mixture_specific_heat_mid_[j+1] -
	       params->thermal_conductivity_[j]*inv_dz[jext]/
	       params->mixture_specific_heat_mid_[j])*
	      params->rel_vol_ext_[jext]*inv_dzm[jext] -
	      c*y_ptr[j*num_states+num_species]*params->rel_vol_ext_[jext];
	  }
	} //for k < num_nonzeros
      }//if fixed T point

      if(params->pseudo_unsteady_) {
        // Add -1/dt term to Yi and T
        for(int k=0; k<num_species; ++k) {
          params->saved_jacobian_chem_[j*num_nonzeros_zerod+params->diagonal_id_chem_[k]] -= 1.0/params->dt_;
        }
        params->saved_jacobian_chem_[j*num_nonzeros_zerod+params->diagonal_id_chem_[num_species+1]] -= 1.0/params->dt_;
      }

      //Add/subtract identity
      for(int k=0; k<num_states; ++k) {
	params->saved_jacobian_chem_[j*num_nonzeros_zerod+params->diagonal_id_chem_[k]] -= constant;
      }

    } // for j<num_local_points

    for(int j=0; j<num_local_points; ++j) {
      for(int k=0; k<num_nonzeros_zerod; ++k) {
	params->reactor_jacobian_chem_[k] =
	  params->saved_jacobian_chem_[j*num_nonzeros_zerod+k];
      }
      // factor the numerical jacobian
      if(params->sparse_matrix_chem_[j]->IsFirstFactor()) {
        error_flag =
          params->sparse_matrix_chem_[j]->FactorNewPatternCCS(num_nonzeros_zerod,
							      &params->row_id_chem_[0],
							      &params->column_sum_chem_[0],
							      &params->reactor_jacobian_chem_[0]);
      } else {
        error_flag =
          params->sparse_matrix_chem_[j]->FactorSamePattern(&params->reactor_jacobian_chem_[0]);
      }
      if(error_flag != 0) {
	printf("Sparse matrix error at point %d\n", j);
//	params->logger_->PrintF(
//				"# DEBUG: grid point %d (z = %.18g [m]) reactor produced a\n"
//				"#        sparse matrix error flag = %d\n",
//				j,
//				params->z_[j],
//				error_flag);
	return error_flag;
      }

    } // for j<num_local_points

  } else {
    // recompute and factor the Jacobian, there is no saved data
    for(int j=0; j<num_local_points; ++j) {
      int jglobal = j + my_pe*num_local_points;
      if(jglobal == params->j_fix_) {
	Tfix = true;
      } else {
	Tfix = false;
      }
      params->reactor_->GetJacobianSteady(&y_ptr[j*num_states],
					  &params->rhsConv_[j*num_states],
					  Tfix,
                                          &params->step_limiter_[0],
					  &params->reactor_jacobian_chem_[0]);

      // dmdot_j/dT_j at j=j_fix
      //int jglobal = j + my_pe*num_local_points;
      int jext = j + nover;
      if(jglobal == params->j_fix_) {
	for(int k=0; k<num_nonzeros_zerod; ++k) {
	  int col_id = params->column_id_chem_[k];
	  int row_id = params->row_id_chem_[k];
	  if (col_id==num_species+1 && row_id==num_species) {
	    double c = (dz[jext+1]-dz[jext])/dz[jext+1]/dz[jext];

	    params->reactor_jacobian_chem_[k] =
	      params->reactor_jacobian_chem_[params->diagonal_id_chem_[num_states-1]] +
	      (-params->thermal_conductivity_[j+1]*inv_dz[jext+1]/
	       params->mixture_specific_heat_mid_[j+1] -
	       params->thermal_conductivity_[j]*inv_dz[jext]/
	       params->mixture_specific_heat_mid_[j])*
	      params->rel_vol_ext_[jext]*inv_dzm[jext] -
	      c*y_ptr[j*num_states+num_species]*params->rel_vol_ext_[jext];
	  }
        }
      }//if fixed T point

      if(params->pseudo_unsteady_) {
        // Add -1/dt term to Yi and T
        for(int k=0; k<num_species; ++k) {
          params->saved_jacobian_chem_[j*num_nonzeros_zerod+params->diagonal_id_chem_[k]] -= 1.0/params->dt_;
        }
        params->saved_jacobian_chem_[j*num_nonzeros_zerod+params->diagonal_id_chem_[num_species+1]] -= 1.0/params->dt_;
      }

      //Add/subtract identity
      for(int k=0; k<num_states; ++k) {
	params->reactor_jacobian_chem_[params->diagonal_id_chem_[k]] -= constant;
      }

      // factor the numerical jacobian
      if(params->sparse_matrix_chem_[j]->IsFirstFactor()) {
        error_flag =
          params->sparse_matrix_chem_[j]->FactorNewPatternCCS(num_nonzeros_zerod,
							      &params->row_id_chem_[0],
							      &params->column_sum_chem_[0],
							      &params->reactor_jacobian_chem_[0]);
      } else {
        error_flag =
          params->sparse_matrix_chem_[j]->FactorSamePattern(&params->reactor_jacobian_chem_[0]);
      }
      if(error_flag != 0) {
	printf("Sparse matrix error flag = %d\n", error_flag);
//	params->logger_->PrintF(
//				"# DEBUG: grid point %d (z = %.18g [m]) reactor produced a\n"
//				"#        sparse matrix error flag = %d\n",
//				j,
//				params->z_[j],
//				error_flag);
	return error_flag;
      }
    } // for(int j=0; j<num_local_points; ++j)
  } // if(params->store_jacobian_)

      // Add/Subtract identity to/from transport jacobian
  for(int j=0; j<num_local_points; ++j) {
      for(int k=0; k<num_states; ++k) {
        params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] += constant;
      }
  }

      // Multiply by inverse of chemical jacobian
      double inverse_chem_jacobian;
  for(int j=0; j<num_local_points; ++j) {
      for(int k=0; k<num_states; ++k) {
      inverse_chem_jacobian = 1.0/params->saved_jacobian_chem_[j*num_nonzeros_zerod+params->diagonal_id_chem_[k]];
        params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] *= inverse_chem_jacobian;
        params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 - 1] *= inverse_chem_jacobian;
        params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 1] *= inverse_chem_jacobian;
      }
  }

  // Add identity
  for(int j=0; j<num_local_points; ++j) {
    for(int k=0; k<num_states; ++k) {
      params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] += 1.0;
    }
  }

#ifdef ZERORK_MPI
  // Communications to solve banded transport Jacobian
  // Each processor handles the full grid for a subset of species
  MPI_Comm comm = params->comm_;
  long int dsize = num_local_points*5;
  int nodeDest;

  // Gather the Jacobian to have all grid points for each species
  for(int j=0; j<num_states; ++j) {
    nodeDest = j/params->num_states_per_proc_;
    int jlocal = j % params->num_states_per_proc_;
    int start_band = j*(num_local_points*5);
    int start_band2 = jlocal*(num_total_points*5);
    double* dest = params->my_pe_ == nodeDest ? &params->banded_jacobian2_[start_band2] : nullptr;

    MPI_Gather(&params->banded_jacobian_[start_band],
    	       dsize,
    	       PVEC_REAL_MPI_TYPE,
               dest,
    	       dsize,
    	       PVEC_REAL_MPI_TYPE,
	       nodeDest,
    	       comm);
  }
#else
  for(int j=0; j<params->banded_jacobian2_.size(); j++) {
    params->banded_jacobian2_[j] = params->banded_jacobian_[j];
  }
#endif
  // Reorder
  for(int j=0; j<num_states_local; ++j) {
    for(int i=0; i<num_total_points; ++i) {
      for(int s=0; s<4; ++s) {
	params->banded_jacobian_serial_[j*(num_total_points*4) + i*4 + s] =
	  params->banded_jacobian2_[j*(num_total_points*5) + i*5 + s + 1];
      }
    }
  }

  // Factorize for each species
  int one = 1;
  int LDAB = 4;
  int dim = num_total_points;
  for(int j=0; j<num_states_local; ++j) {
    dgbtrf_(&dim,
	    &dim,
	    &one,
	    &one,
	    &params->banded_jacobian_serial_[j*num_total_points*4],
	    &LDAB,
	    &params->pivots_serial_[j*num_total_points],
	    &error_flag);
  }

  return error_flag;
}

#if defined SUNDIALS2
int ReactorAFSolve(N_Vector y, // [in] state vector
		   N_Vector yscale, // [in] state scaler
		   N_Vector ydot, // [in] state derivative
		   N_Vector ydotscale, // [in] state derivative scaler
		   N_Vector vv, // [in/out] rhs/solution vector
		   void *user_data, // [in/out]
                   N_Vector tmp)
{
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorAFSolve(N_Vector y, // [in] state vector
		   N_Vector yscale, // [in] state scaler
		   N_Vector ydot, // [in] state derivative
		   N_Vector ydotscale, // [in] state derivative scaler
		   N_Vector vv, // [in/out] rhs/solution vector
		   void *user_data) // [in/out]
{
#endif
  double *solution = N_VGetArrayPointer(vv);
  int error_flag = 0;

  error_flag = AFSolve(&solution[0], user_data);
  return error_flag;

}

// Solve approximately factorized Jacobian
int AFSolve(double solution[],
            void *user_data)
{
  FlameParams *params = (FlameParams *)user_data;
  const int num_total_points  = params->num_points_;
  const int num_local_points  = params->num_local_points_;
  const int num_states  = params->reactor_->GetNumStates();
  const int num_states_local = params->num_states_local_;
  int error_flag = 0;

  // Solve Local sparse chemistry with SuperLU
  int start_id=0;
  for(int j=0; j<num_local_points; ++j) {
    error_flag = params->sparse_matrix_chem_[j]->Solve(&solution[start_id],
						       &solution[start_id]);
    start_id += num_states;
    if(error_flag != 0) {
      printf("AFSolve sparse matrix error: %d\n", error_flag);
      return error_flag;
    }
  }

  // Banded transport
  // Communications for banded_jacobian2
#ifdef ZERORK_MPI
  MPI_Comm comm = params->comm_;
#endif
  long int dsize = num_local_points;
  int nodeDest, nodeFrom;

  std::vector<double> solution_allspecies, solution_species;
  solution_allspecies.assign(num_total_points*num_states_local, 0.0);
  solution_species.assign(num_local_points*num_states, 0.0);

  // Reorder solution vector by species
  for(int j=0; j<num_states; ++j)
    for(int i=0; i<num_local_points; ++i)
      solution_species[j*num_local_points+i] = solution[j+i*num_states];

#ifdef ZERORK_MPI
  // Gather all grid points for each species
  for(int j=0; j<num_states; ++j) {
    nodeDest = j/params->num_states_per_proc_;
    int jlocal = j % params->num_states_per_proc_;
    int start_id = j*num_local_points;
    int start_id2 = jlocal*num_total_points;
    double* dest = params->my_pe_ == nodeDest ? &solution_allspecies[start_id2] : nullptr;

    MPI_Gather(&solution_species[start_id],
    	       dsize,
    	       PVEC_REAL_MPI_TYPE,
               dest,
    	       dsize,
    	       PVEC_REAL_MPI_TYPE,
	       nodeDest,
    	       comm);
  }
#else
  for(int j=0; j<solution_allspecies.size(); j++) {
    solution_allspecies[j] = solution_species[j];
  }
#endif
  // Solve banded matrix for each species
  int dim = num_total_points;
  int one = 1;
  int LDAB = 4;
  int LDB = num_total_points;
  for(int j=0; j<num_states_local; ++j) {
    dgbtrs_("N",
	    &dim,
	    &one,
	    &one,
	    &one,
	    &params->banded_jacobian_serial_[j*num_total_points*4],
	    &LDAB,
	    &params->pivots_serial_[j*num_total_points],
	    &solution_allspecies[j*num_total_points],
	    &LDB,
	    &error_flag);

    if(error_flag != 0)
      printf("AFSolve banded matrix error: %d\n", error_flag);
  }
#ifdef ZERORK_MPI
  //Scatter back the solution vector for each species
  for(int j=0; j<num_states; ++j) {
    nodeFrom = j/params->num_states_per_proc_;
    int jlocal = j % params->num_states_per_proc_;
    int start_id = j*num_local_points;
    int start_id2 = jlocal*num_total_points;
    double* source = params->my_pe_ == nodeFrom ? &solution_allspecies[start_id2] : nullptr;

    MPI_Scatter(source,
    		dsize,
    		PVEC_REAL_MPI_TYPE,
    		&solution_species[start_id],
    		dsize,
    		PVEC_REAL_MPI_TYPE,
    		nodeFrom,
    		comm);
  }
#else
  for(int j=0; j<solution_allspecies.size(); j++)
    solution_species[j] = solution_allspecies[j];
#endif
  // Reorder solution vector by grid points
  for(int j=0; j<num_states; ++j)
    for(int i=0; i<num_local_points; ++i)
      solution[j+i*num_states] = solution_species[j*num_local_points+i];

  return error_flag;
}


