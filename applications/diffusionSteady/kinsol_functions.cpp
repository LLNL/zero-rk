#include "kinsol_functions.h"
#include "flame_params.h"
#include "utilities/math_utilities.h"

extern "C" void dgbtrf_(int* dim1, int* dim2, int* nu, int* nl, double* a, int* lda, int* ipiv, int* info);
extern "C" void dgbtrs_(char *TRANS, int *N, int *NRHS, int* nu, int* nl, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

// Upwind scheme for convective term
static double NonLinearConvectUpwind(double velocity,
                                     double y_previous,
                                     double y_current,
                                     double y_next,
                                     double inv_dz_prev,
                                     double inv_dz);

// Main RHS function
int ConstPressureFlame(N_Vector y,
		       N_Vector ydot,
		       void *user_data)
{
  FlameParams *params = (FlameParams *)user_data;
  const int num_local_points = params->num_local_points_;
  const int num_states  = params->reactor_->GetNumStates();
  int Nlocal = num_local_points*num_states;

  ConstPressureFlameComm(Nlocal, y, user_data);
  ConstPressureFlameLocal(Nlocal, y, ydot, user_data);

  return 0;
}

// Parallel communications using MPI
int ConstPressureFlameComm(int nlocal,
                           N_Vector y,
                           void *user_data)
{
  FlameParams *params = (FlameParams *)user_data;
  double *y_ptr    = NV_DATA_P(y);
  const int num_local_points = params->num_local_points_;
  const int num_states  = params->reactor_->GetNumStates();
  int my_pe = params->my_pe_;
  int npes  = params->npes_;

  MPI_Comm comm = params->comm_;
  MPI_Status status;
  const int nover = params->nover_;
  long int dsize = num_states*nover;

  // Copy y_ptr data into larger arrays
  for (int j=0; j<num_states*num_local_points; ++j)
    params->y_ext_[num_states*nover + j] = y_ptr[j];

  // Send/receive
  int nodeDest = my_pe-1;
  if (nodeDest < 0) nodeDest = npes-1;
  int nodeFrom = my_pe+1;
  if (nodeFrom > npes-1) nodeFrom = 0;
  MPI_Sendrecv(&params->y_ext_[nover*num_states],
               dsize, PVEC_REAL_MPI_TYPE,
               nodeDest, 0,
               &params->y_ext_[num_states*(num_local_points+nover)],
               dsize, PVEC_REAL_MPI_TYPE,
               nodeFrom, 0, comm, &status);

  nodeDest = my_pe+1;
  if (nodeDest > npes-1) nodeDest = 0;
  nodeFrom = my_pe-1;
  if (nodeFrom < 0) nodeFrom = npes-1;
  MPI_Sendrecv(&params->y_ext_[num_states*num_local_points],
               dsize, PVEC_REAL_MPI_TYPE,
               nodeDest, 0,
               &params->y_ext_[0],
               dsize, PVEC_REAL_MPI_TYPE,
               nodeFrom, 0, comm, &status);

  return 0;
}

int ConstPressureFlameLocal(int nlocal,
			    N_Vector y,
			    N_Vector ydot,
			    void *user_data)
{
  FlameParams *params = (FlameParams *)user_data;
  double *y_ptr    = NV_DATA_P(y);   // caution: assumes realtype == double
  double *ydot_ptr = NV_DATA_P(ydot); // caution: assumes realtype == double
  const int num_total_points = params->z_.size();
  const int num_local_points = params->num_local_points_;
  const int num_states  = params->reactor_->GetNumStates();
  const int num_species = params->reactor_->GetNumSpecies();
  const int num_local_states = num_local_points*num_states;
  const int nover = params->nover_;
  int my_pe = params->my_pe_;
  int npes  = params->npes_;
  MPI_Comm comm = params->comm_;

  const double ref_temperature = params->ref_temperature_;

  std::vector<double> enthalpies;
  enthalpies.assign(num_species,0.0);

  std::vector<double> rhs;
  rhs.assign(num_local_points*num_states,0.0);

  std::vector<double> conductivity_over_cp,
    dissipation_rate_times_rho, mass_fraction_over_mixture_mass,
    mass_fraction_over_mixture_mass_times_sum;
  dissipation_rate_times_rho.assign(num_local_points+2*nover, 0.0);
  conductivity_over_cp.assign(num_local_points+2*nover, 0.0);
  mass_fraction_over_mixture_mass.assign((num_local_points+2*nover)*num_species, 0.0);
  mass_fraction_over_mixture_mass_times_sum.assign((num_local_points+2*nover)*num_species, 0.0);

  double relative_volume_j, convection_velocity_j;
  int transport_error;

  double local_max;
  double thermal_diffusivity;

  // Set the RHS to zero
  for(int j=0; j<num_local_states; ++j) {
    ydot_ptr[j] = 0.0;
    rhs[j] = 0.0;
  }

  // compute the constant pressure reactor source term
  if(params->fix_temperature_) {
    // Different function if temperature is fixed
    for(int j=0; j<num_local_points; ++j)
      params->reactor_->GetTimeDerivativeDiffusionSteadyFixT(
        &y_ptr[j*num_states],
        &params->step_limiter_[0],
        &rhs[j*num_states]);
  } else {
    for(int j=0; j<num_local_points; ++j)
      params->reactor_->GetTimeDerivativeDiffusionSteady(
        &y_ptr[j*num_states],
        &params->step_limiter_[0],
        &rhs[j*num_states]);
  }

  // compute soot
  if(params->soot_) {
    for(int j=0; j<num_local_points; ++j) {
      int jext = j + nover;

      CalcRhoDot(params,
                 &y_ptr[j*num_states],
                 params->rho_dot_[jext]);

      UpdateProductionRates(params,
                            &y_ptr[j*num_states],
                            &rhs[j*num_states]);
    }
  }

  //--------------------------------------------------------------------------
  std::vector<double> dz, inv_dz;
  dz.assign( num_local_points+(2*nover), 0.0);
  inv_dz.assign( num_local_points+(2*nover), 0.0);

  for (int j=0; j<num_local_points+2*nover; ++j) {
    dz[j] = params->dz_local_[j];
    inv_dz[j] = params->inv_dz_local_[j];
  }

  //--------------------------------------------------------------------------
  // Boundary Conditions
  // First proc: oxidizer conditions in ghost cells (Z=0)
  if (my_pe ==0) {
    for(int j=0; j<nover; ++j) {
      for(int k=0; k<num_species; ++k) {
	params->y_ext_[j*num_states + k] = params->oxidizer_mass_fractions_[k];
      }
      params->y_ext_[(j+1)*num_states-1] = params->oxidizer_temperature_;
      params->y_ext_[(j+1)*num_states-2] = params->oxidizer_relative_volume_;
    }
  }

  // Last proc: fuel conditions in ghost cells (Z=1)
  if (my_pe == npes-1) {
    for(int j=num_local_points+nover; j<num_local_points+2*nover; ++j) {
      for(int k=0; k<num_species; ++k) {
	params->y_ext_[j*num_states + k] = params->fuel_mass_fractions_[k];
      }
      params->y_ext_[(j+1)*num_states-1] = params->fuel_temperature_;
      params->y_ext_[(j+1)*num_states-2] = params->fuel_relative_volume_;
    }
  }

  //--------------------------------------------------------------------------
  // Compute the interior heat capacity, conductivity, species mass fluxes/Lewis,
  // mixture molecular weight, dissipation rate, etc.
  // including ghost cells
  for(int j=0; j<num_local_points+2*nover; ++j) { //+2
    int jlocal = j-nover;
    int jext = j;
    int jglobal = jlocal + my_pe*num_local_points;

    relative_volume_j = params->y_ext_[jext*num_states+num_species];

    // mass fractions
    for(int k=0; k<num_species; ++k)
      params->transport_input_.mass_fraction_[k] = params->y_ext_[jext*num_states+k];

    // temperature
    params->transport_input_.temperature_ = ref_temperature*
      params->y_ext_[(jext+1)*num_states-1];

    // specific heat at grid point j
    params->mixture_specific_heat_[j] =
      params->reactor_->GetMixtureSpecificHeat_Cp(
        ref_temperature*params->y_ext_[(jext+1)*num_states-1],
        &params->y_ext_[jext*num_states],
        &params->species_specific_heats_[num_species*j]);

    // compute the conductivity at grid point j
    transport_error = params->trans_->GetMixtureConductivity(
      params->transport_input_,
      &params->thermal_conductivity_[j]);
    if(transport_error != transport::NO_ERROR) {
      return transport_error;
    }

    // compute the species mass flux at the upstream mid point
    // only used to get Lewis numbers here
    // should be a different function
    transport_error = params->trans_->GetSpeciesMassFlux(
      params->transport_input_,
      num_species,
      nullptr,
      nullptr,
      &params->species_mass_flux_[j*num_species],
      &params->species_lewis_numbers_[j*num_species]);
    if(transport_error != transport::NO_ERROR) {
      return transport_error;
    }

    // compute mixture molecular weight
    double mass_fraction_weight_sum = 0.0;
    for(int k=0; k<num_species; ++k) {
      mass_fraction_weight_sum += params->inv_molecular_mass_[k]*
        params->y_ext_[jext*num_states+k];
    }
    params->mixture_molecular_mass_[j] = 1.0/mass_fraction_weight_sum;

    // compute mixture conductivity over heat capacitiy
    conductivity_over_cp[j] = params->thermal_conductivity_[j]/
      params->mixture_specific_heat_[j];

    // compute dissipation rate times density
    dissipation_rate_times_rho[j] = params->dissipation_rate_[j]/relative_volume_j;

    // compute sum of mass fraction over Lewis number
    double local_sum_Y_over_Le = 0.0;
    for(int k=0; k<num_species; ++k) {
      local_sum_Y_over_Le += params->y_ext_[jext*num_states+k]/
        params->species_lewis_numbers_[j*num_species + k];
    }
    params->sum_mass_fraction_over_Lewis_[j] = local_sum_Y_over_Le;

    // compute mass fraction over mixture molecular mass and
    // mass fraction over mixture mass times sum of mass fraction over Lewis
    for(int k=0; k<num_species; ++k) {
      mass_fraction_over_mixture_mass[j*num_species+k] =
        params->y_ext_[jext*num_states+k]/
        params->mixture_molecular_mass_[j];
      mass_fraction_over_mixture_mass_times_sum[j*num_species+k] =
        mass_fraction_over_mixture_mass[j*num_species+k]*
        params->sum_mass_fraction_over_Lewis_[j];
    }

  } // for j<num_local_points+2*nover


  //--------------------------------------------------------------------------
  // Pre-compute derivatives of sum(Y_i/Le_i)
  for(int j=0; j<num_local_points; ++j) {
    int jext = j + nover;

    // First derivative
    //coefficients of j+1, j, j-1 terms
    double b=0,c=0,d=0;
    b = dz[jext]/dz[jext+1]/(dz[jext]+dz[jext+1]);
    c = (dz[jext+1]-dz[jext])/dz[jext+1]/dz[jext];
    d = -dz[jext+1]/dz[jext]/(dz[jext]+dz[jext+1]);

    // Second derivative
    // Centered three-point stencil for now
    double bb=0, cc=0, dd=0;
    double denom = dz[jext]*dz[jext+1]*(dz[jext]+dz[jext+1]);
    bb =  2.0*dz[jext]/denom;
    cc = -2.0*(dz[jext]+dz[jext+1])/denom;
    dd =  2.0*dz[jext+1]/denom;

    params->sum_mass_fraction_over_Lewis_grad_[jext] = 0.0;
    params->sum_mass_fraction_over_Lewis_laplacian_[jext] = 0.0;

    for(int k=0; k<num_species; ++k){
      params->sum_mass_fraction_over_Lewis_grad_[jext] +=
        (b*params->y_ext_[(jext+1)*num_states+k] +
         c*params->y_ext_[jext*num_states+k] +
         d*params->y_ext_[(jext-1)*num_states+k])/
        params->species_lewis_numbers_[jext*num_species + k];

      params->sum_mass_fraction_over_Lewis_laplacian_[jext] +=
        (bb*params->y_ext_[(jext+1)*num_states+k] +
         cc*params->y_ext_[jext*num_states+k] +
         dd*params->y_ext_[(jext-1)*num_states+k])/
        params->species_lewis_numbers_[jext*num_species + k];
    }
  }

  //--------------------------------------------------------------------------
  // Compute convective and diffusive terms for species and temperature
  for(int j=0; j<num_local_points; ++j) {
    int jext = j + nover;
    int jglobal = j + my_pe*num_local_points;

    // First derivative
    // Only centered for now
    //coefficients of j+2, j+1, j, j-1, j-2 terms
    double a=0,b=0,c=0,d=0,e=0;
    a = 0;
    b = dz[jext]/dz[jext+1]/(dz[jext]+dz[jext+1]);
    c = (dz[jext+1]-dz[jext])/dz[jext+1]/dz[jext];
    d = -dz[jext+1]/dz[jext]/(dz[jext]+dz[jext+1]);
    e = 0;

    // Second derivative
    // Centered three-point stencil for now
    double bb=0, cc=0, dd=0;
    double denom = dz[jext]*dz[jext+1]*(dz[jext]+dz[jext+1]);
    bb =  2.0*dz[jext]/denom;
    cc = -2.0*(dz[jext]+dz[jext+1])/denom;
    dd =  2.0*dz[jext+1]/denom;

    relative_volume_j = params->y_ext_[jext*num_states+num_species];

    // compute the convection velocity
    if(params->fix_temperature_) {
      params->convection_velocity_[jext] = -270.0*pow(params->z_[jglobal],1.8)*
        pow(1.0-params->z_[jglobal],0.75);
    } else {
      params->convection_velocity_[jext] = 0.25*relative_volume_j*(
        (b*dissipation_rate_times_rho[jext+1] +
         c*dissipation_rate_times_rho[jext] +
         d*dissipation_rate_times_rho[jext-1]) +
        dissipation_rate_times_rho[jext]/conductivity_over_cp[jext]*
        (b*conductivity_over_cp[jext+1] +
         c*conductivity_over_cp[jext] +
         d*conductivity_over_cp[jext-1]));
    }

    convection_velocity_j = params->convection_velocity_[jext];

    // compute the species convection/diffusion terms
    for(int k=0; k<num_species; ++k) {

      // Diffusion - 2: 1/2*chi/Le * ddY/ddZ
      rhs[j*num_states+k] += 0.5*params->dissipation_rate_[jext]/
        params->species_lewis_numbers_[jext*num_species + k]*
        (bb*params->y_ext_[(jext+1)*num_states+k] +
         cc*params->y_ext_[jext*num_states+k] +
         dd*params->y_ext_[(jext-1)*num_states+k]);

      // Non-unity Lewis terms below
      if(!params->unity_Lewis_) {
        // Following three terms are hardest to converge
        if(params->full_equations_) {

          // Diffusion - 4: 1/2*chi/Le*ddW/ddZ
          params->mixture_molecular_mass_laplacian_[jext] =
            (bb*params->mixture_molecular_mass_[jext+1] +
             cc*params->mixture_molecular_mass_[jext] +
             dd*params->mixture_molecular_mass_[jext-1]);

          rhs[j*num_states+k] += 0.5*params->dissipation_rate_[jext]/
            params->species_lewis_numbers_[jext*num_species + k]*
            params->y_ext_[jext*num_states+k]/
            params->mixture_molecular_mass_[jext]*
            params->mixture_molecular_mass_laplacian_[jext];

          // Diffusion correction - 5.1: 1/2*chi*X*dd(sumYi/Lei)/ddZ
          rhs[j*num_states+k] -= 0.5*params->dissipation_rate_[jext]*
            params->y_ext_[jext*num_states+k]*
            params->sum_mass_fraction_over_Lewis_laplacian_[jext];

          // Diffusion correction - 5.2: 1/2*chi*Yi/W*sumYi/Lei*ddW/ddZ
          rhs[j*num_states+k] -= 0.5*params->dissipation_rate_[jext]*
            params->y_ext_[jext*num_states+k]/
            params->mixture_molecular_mass_[jext]*
            params->sum_mass_fraction_over_Lewis_[jext]*
            params->mixture_molecular_mass_laplacian_[jext];
        }

        // Mass fraction convection velocity - 6: u*(1/Le-1)*dY/dZ
        double uConv = convection_velocity_j*
          (1.0/params->species_lewis_numbers_[jext*num_species + k] - 1.0);

        rhs[j*num_states+k] -= NonLinearConvectUpwind(-uConv,
                                                      params->y_ext_[(jext-1)*num_states+k],
                                                      params->y_ext_[jext*num_states+k],
                                                      params->y_ext_[(jext+1)*num_states+k],
                                                      inv_dz[jext],
                                                      inv_dz[jext+1]);

        // Molar mass convection velocity - 7: dW/dZ*1/Le*(u*Yi/W + chi/2*d(Yi/W)/dZ)
        params->mixture_molecular_mass_grad_[jext] =
          (b*params->mixture_molecular_mass_[jext+1] +
           c*params->mixture_molecular_mass_[jext] +
           d*params->mixture_molecular_mass_[jext-1]);

        rhs[j*num_states+k] += params->mixture_molecular_mass_grad_[jext]/
          params->species_lewis_numbers_[jext*num_species + k]*(
            convection_velocity_j*params->y_ext_[jext*num_states+k]/
            params->mixture_molecular_mass_[jext] +
            0.5*params->dissipation_rate_[jext]*
            (b*mass_fraction_over_mixture_mass[(jext+1)*num_species+k] +
             c*mass_fraction_over_mixture_mass[(jext)*num_species+k] +
             d*mass_fraction_over_mixture_mass[(jext-1)*num_species+k]));

        // Mass fraction convection velocity - 8: d(sumYi/Lei)/dZ*(u*Yi+chi/2*dYi/dZ)
        rhs[j*num_states+k] -= params->sum_mass_fraction_over_Lewis_grad_[jext]*(
          convection_velocity_j*params->y_ext_[jext*num_states+k] +
          0.5*params->dissipation_rate_[jext]*
          (a*params->y_ext_[(jext+2)*num_states+k] +
           b*params->y_ext_[(jext+1)*num_states+k] +
           c*params->y_ext_[jext*num_states+k] +
           d*params->y_ext_[(jext-1)*num_states+k] +
           e*params->y_ext_[(jext-2)*num_states+k]));

        // Molar mass convection correction - 9:
        rhs[j*num_states+k] -= params->mixture_molecular_mass_grad_[jext]*(
          convection_velocity_j*params->y_ext_[jext*num_states+k]/
          params->mixture_molecular_mass_[jext]*
          params->sum_mass_fraction_over_Lewis_[jext] +
          0.5*params->dissipation_rate_[jext]*
          (b*mass_fraction_over_mixture_mass_times_sum[(jext+1)*num_species+k] +
           c*mass_fraction_over_mixture_mass_times_sum[(jext)*num_species+k] +
           d*mass_fraction_over_mixture_mass_times_sum[(jext-1)*num_species+k]));

      } //if(!unity_Lewis)

      if(params->soot_) {
        // -rhodot*(Y_k - Z*dY_k/dZ)
        rhs[j*num_states+k] -= params->rho_dot_[jext]*relative_volume_j*
          (params->y_ext_[jext*num_states+k] - params->z_[jglobal]*
           (a*params->y_ext_[(jext+2)*num_states+k] +
            b*params->y_ext_[(jext+1)*num_states+k] +
            c*params->y_ext_[jext*num_states+k] +
            d*params->y_ext_[(jext-1)*num_states+k] +
            e*params->y_ext_[(jext-2)*num_states+k]));
      }
    } // for k < num_species

    // Temperature equation
    if(params->fix_temperature_) {
      rhs[(j+1)*num_states-1] = params->y_ext_[(jext+1)*num_states-1] -
        params->fixed_temperature_[jglobal];
    } else {
      // Heat conduction - 2
      rhs[(j+1)*num_states-1] += 0.5*params->dissipation_rate_[jext]*
        (bb*params->y_ext_[(jext+2)*num_states-1] +
         cc*params->y_ext_[(jext+1)*num_states-1] +
         dd*params->y_ext_[jext*num_states-1]);

      // Heat condution - 3
      rhs[(j+1)*num_states-1] += 0.5*params->dissipation_rate_[jext]/
        params->mixture_specific_heat_[jext]*
        (b*params->mixture_specific_heat_[jext+1] +
         c*params->mixture_specific_heat_[jext] +
         d*params->mixture_specific_heat_[jext-1])*
        (b*params->y_ext_[(jext+2)*num_states-1] +
         c*params->y_ext_[(jext+1)*num_states-1] +
         d*params->y_ext_[jext*num_states-1]);

      // Enthalpy flux - 7
      double enthalpy_flux_sum = 0.0;
      for(int k=0; k<num_species; k++) {

        enthalpy_flux_sum -= params->species_specific_heats_[num_species*jext + k]/
          params->species_lewis_numbers_[jext*num_species + k]*
          (a*params->y_ext_[(jext+2)*num_states+k] +
           b*params->y_ext_[(jext+1)*num_states+k] +
           c*params->y_ext_[jext*num_states+k] +
           d*params->y_ext_[(jext-1)*num_states+k] +
           e*params->y_ext_[(jext-2)*num_states+k]);

        if(!params->unity_Lewis_) {

          enthalpy_flux_sum += params->mixture_specific_heat_[jext]/
            params->species_lewis_numbers_[jext*num_species + k]*
            (a*params->y_ext_[(jext+2)*num_states+k] +
             b*params->y_ext_[(jext+1)*num_states+k] +
             c*params->y_ext_[jext*num_states+k] +
             d*params->y_ext_[(jext-1)*num_states+k] +
             e*params->y_ext_[(jext-2)*num_states+k]);

          enthalpy_flux_sum += (params->mixture_specific_heat_[jext] -
                                params->species_specific_heats_[jext*num_species + k])/
            params->species_lewis_numbers_[jext*num_species + k]*
            params->y_ext_[jext*num_states+k]/
            params->mixture_molecular_mass_[jext]*
            params->mixture_molecular_mass_grad_[jext];

        }// if !unity_Lewis
      } // for k<num_species

      // Save enthalpy_flux_sum for use in jacobian
      params->enthalpy_flux_sum_[jext] = enthalpy_flux_sum;

      rhs[(j+1)*num_states-1] -= 0.5*params->dissipation_rate_[jext]/
        params->mixture_specific_heat_[jext]*
        params->enthalpy_flux_sum_[jext]*
        (b*params->y_ext_[(jext+2)*num_states-1] +
         c*params->y_ext_[(jext+1)*num_states-1] +
         d*params->y_ext_[jext*num_states-1]);

      if(params->soot_) {
        // rhodot*Z*dcpT/dZ (divided by rho*cp for correct units)
        rhs[(j+1)*num_states-1] += params->rho_dot_[jext]*params->z_[jglobal]*relative_volume_j/
          params->mixture_specific_heat_[jext]*
          (b*params->y_ext_[(jext+2)*num_states-1]*params->mixture_specific_heat_[jext+1] +
           c*params->y_ext_[(jext+1)*num_states-1]*params->mixture_specific_heat_[jext] +
           d*params->y_ext_[jext*num_states-1]*params->mixture_specific_heat_[jext-1]);
      }

    } // if fix_temperature
  } // for(int j=0; j<num_local_points; ++j) // loop computing rhs

  // -------------------------------------------------------------------------
  // Compute the final residual
  for(int j=0; j<num_local_points; ++j) {
    int rvol_id = j*num_states+num_species; // relative volume index of pt j
    int temp_id = (j+1)*num_states-1;  // temperature index of pt j

    for(int k=0; k<num_species; ++k) {
      ydot_ptr[j*num_states+k] = rhs[j*num_states+k];
    }
    ydot_ptr[rvol_id] = rhs[rvol_id]; //handled in GetTimeDerivative
    ydot_ptr[temp_id] = rhs[temp_id];
  }

  // Add time derivative term if pseudo unsteady
  if(params->pseudo_unsteady_) {
    for(int j=0; j<num_local_points; ++j) {
      int rvol_id = j*num_states+num_species; // relative volume index of pt j
      int temp_id = rvol_id+1 ;               // temperature index of pt j
      int mom_id  = rvol_id+2;                // momentum index of pt j

      for(int k=0; k<num_species; ++k) {
        ydot_ptr[j*num_states+k] -= (y_ptr[j*num_states+k] - params->y_old_[j*num_states+k])/
          params->dt_;
      }
      ydot_ptr[temp_id] -= (y_ptr[temp_id] - params->y_old_[temp_id])/params->dt_;
    }
  }

  // -------------------------------------------------------------------------
  // Parallel communication for finite difference Jacobian
  // Move to a separate function?
  if(params->integrator_type_ == 2) {
    MPI_Status status;
    long int dsize = num_states*nover;

    // Copy ydot_ptr into larger arrays
    for (int j=0; j<num_states*num_local_points; ++j)
      params->rhs_ext_[num_states*nover + j] = ydot_ptr[j];

    if(npes>1) {
      // MPI sendrecv
      int nodeDest = my_pe-1;
      if (nodeDest < 0) nodeDest = npes-1;
      int nodeFrom = my_pe+1;
      if (nodeFrom > npes-1) nodeFrom = 0;
      MPI_Sendrecv(&params->rhs_ext_[nover*num_states],
                   dsize, PVEC_REAL_MPI_TYPE,
                   nodeDest, 0,
                   &params->rhs_ext_[num_states*(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE,
                   nodeFrom, 0, comm, &status);

      nodeDest = my_pe+1;
      if (nodeDest > npes-1) nodeDest = 0;
      nodeFrom = my_pe-1;
      if (nodeFrom < 0) nodeFrom = npes-1;
      MPI_Sendrecv(&params->rhs_ext_[num_states*num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE,
                   nodeDest, 0,
                   &params->rhs_ext_[0],
                   dsize, PVEC_REAL_MPI_TYPE,
                   nodeFrom, 0, comm, &status);
    }
  }

  // Parallel communication for analytical transport Jacobian
  if(params->integrator_type_ == 3) {
    MPI_Status status;
    long int dsize = nover;

    if(npes>1) {
      // MPI sendrecv
      int nodeDest = my_pe-1;
      if (nodeDest < 0) nodeDest = npes-1;
      int nodeFrom = my_pe+1;
      if (nodeFrom > npes-1) nodeFrom = 0;
      MPI_Sendrecv(&params->enthalpy_flux_sum_[nover],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->enthalpy_flux_sum_[(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->mixture_molecular_mass_[nover],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->mixture_molecular_mass_[(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->mixture_molecular_mass_[nover],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->mixture_molecular_mass_[(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->mixture_molecular_mass_grad_[nover],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->mixture_molecular_mass_grad_[(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->mixture_molecular_mass_laplacian_[nover],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->mixture_molecular_mass_laplacian_[(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->sum_mass_fraction_over_Lewis_[nover],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->sum_mass_fraction_over_Lewis_[(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->sum_mass_fraction_over_Lewis_[nover],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->sum_mass_fraction_over_Lewis_[(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->sum_mass_fraction_over_Lewis_grad_[nover],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->sum_mass_fraction_over_Lewis_grad_[(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->sum_mass_fraction_over_Lewis_laplacian_[nover],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->sum_mass_fraction_over_Lewis_laplacian_[(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->convection_velocity_[nover],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->convection_velocity_[(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->rho_dot_[nover],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->rho_dot_[(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);

      nodeDest = my_pe+1;
      if (nodeDest > npes-1) nodeDest = 0;
      nodeFrom = my_pe-1;
      if (nodeFrom < 0) nodeFrom = npes-1;
      MPI_Sendrecv(&params->enthalpy_flux_sum_[num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->enthalpy_flux_sum_[0],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->mixture_molecular_mass_[num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->mixture_molecular_mass_[0],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->mixture_molecular_mass_[num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->mixture_molecular_mass_[0],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->mixture_molecular_mass_grad_[num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->mixture_molecular_mass_grad_[0],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->mixture_molecular_mass_laplacian_[num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->mixture_molecular_mass_laplacian_[0],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->sum_mass_fraction_over_Lewis_[num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->sum_mass_fraction_over_Lewis_[0],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->sum_mass_fraction_over_Lewis_[num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->sum_mass_fraction_over_Lewis_[0],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->sum_mass_fraction_over_Lewis_grad_[num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->sum_mass_fraction_over_Lewis_grad_[0],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->sum_mass_fraction_over_Lewis_laplacian_[num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->sum_mass_fraction_over_Lewis_laplacian_[0],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->convection_velocity_[num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->convection_velocity_[0],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
      MPI_Sendrecv(&params->rho_dot_[num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
                   &params->rho_dot_[0],
                   dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);
    }// if(npes>1)

    double mixture_molecular_mass_oxidizer=0.0;
    double sum_mass_fraction_over_Lewis_oxidizer=0.0;
    for(int k=0; k<num_species; k++) {
      mixture_molecular_mass_oxidizer += params->oxidizer_mass_fractions_[k]*
        params->inv_molecular_mass_[k];
      sum_mass_fraction_over_Lewis_oxidizer += params->oxidizer_mass_fractions_[k]/
        params->species_lewis_numbers_[k];
    }
    mixture_molecular_mass_oxidizer = 1.0/mixture_molecular_mass_oxidizer;

    double mixture_molecular_mass_fuel=0.0;
    double sum_mass_fraction_over_Lewis_fuel=0.0;
    for(int k=0; k<num_species; k++) {
      mixture_molecular_mass_fuel += params->fuel_mass_fractions_[k]*
        params->inv_molecular_mass_[k];
      sum_mass_fraction_over_Lewis_fuel += params->fuel_mass_fractions_[k]/
        params->species_lewis_numbers_[k];
    }
    mixture_molecular_mass_fuel = 1.0/mixture_molecular_mass_fuel;

    // Boundary Conditions
    // First proc: oxidizer conditions in ghost cells (Z=0)
    if (my_pe ==0) {
      for(int j=0; j<nover; ++j) {
        params->enthalpy_flux_sum_[j] = 0.0;
        params->mixture_molecular_mass_[j] =
          mixture_molecular_mass_oxidizer;
        params->mixture_molecular_mass_[j] =
          mixture_molecular_mass_oxidizer;
        params->mixture_molecular_mass_grad_[j] = 0.0;
        params->mixture_molecular_mass_laplacian_[j] = 0.0;
        params->sum_mass_fraction_over_Lewis_[j] =
          sum_mass_fraction_over_Lewis_oxidizer;
        params->sum_mass_fraction_over_Lewis_[j] =
          sum_mass_fraction_over_Lewis_oxidizer;
        params->sum_mass_fraction_over_Lewis_grad_[j] = 0.0;
        params->sum_mass_fraction_over_Lewis_laplacian_[j] = 0.0;
        params->convection_velocity_[j] = params->convection_velocity_[nover];
        params->rho_dot_[j] = 0.0;
      }
    }

    // Last proc: fuel conditions in ghost cells (Z=1)
    if (my_pe == npes-1) {
      for(int j=num_local_points+nover; j<num_local_points+2*nover; ++j) {
        params->enthalpy_flux_sum_[j] = 0.0;
        params->mixture_molecular_mass_[j] =
          mixture_molecular_mass_fuel;
        params->mixture_molecular_mass_[j] =
          mixture_molecular_mass_fuel;
	params->mixture_molecular_mass_grad_[j] = 0.0;
        params->mixture_molecular_mass_laplacian_[j] = 0.0;
	params->sum_mass_fraction_over_Lewis_[j] =
          sum_mass_fraction_over_Lewis_fuel;
        params->sum_mass_fraction_over_Lewis_[j] =
          sum_mass_fraction_over_Lewis_fuel;
        params->sum_mass_fraction_over_Lewis_grad_[j] = 0.0;
        params->sum_mass_fraction_over_Lewis_laplacian_[j] = 0.0;
	params->convection_velocity_[j] =
          params->convection_velocity_[num_local_points+nover-1];
        params->rho_dot_[j] = 0.0;
      }
    }

  } // if integrator type == 3

  // -------------------------------------------------------------------------
  // initialize the max variables used for explicit time step information
  // compute the max velocity from the mass flux and relative volume stored
  // in the state vector
  double local_temperature;
  local_max = 0.0;
  for(int j=0; j<num_local_points; ++j) {
    int jext = j + nover;
    local_temperature = ref_temperature*params->y_ext_[(jext+1)*num_states-1];
    if(local_temperature > local_max) {
      local_max = local_temperature;
    }
  }
  MPI_Allreduce(&local_max,&params->max_temperature_,1,PVEC_REAL_MPI_TYPE,MPI_MAX,comm);

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
  MPI_Allreduce(&local_max,&params->flame_thickness_,1,PVEC_REAL_MPI_TYPE,MPI_MAX,comm);

  // compute the max thermal diffusivity using the average value of the
  // conductivity and the up and downstream interfaces
  local_max = 0.0;
  for(int j=0; j<num_local_points; ++j) {
    thermal_diffusivity =
      fabs(0.5*(params->thermal_conductivity_[j]+
                params->thermal_conductivity_[j+1])*
           y_ptr[(j+1)*num_states-2]/params->mixture_specific_heat_[j]);
    if(thermal_diffusivity > local_max) {
      local_max = thermal_diffusivity;
    }
  }
  MPI_Allreduce(&local_max,&params->max_thermal_diffusivity_,1,PVEC_REAL_MPI_TYPE,MPI_MAX,comm);

  return 0;
}

// -------------------------------------------------------------------------
// Approximate factorization preconditioner
#if defined SUNDIALS2
int ReactorAFSetup(N_Vector y,      // [in]  state vector
                   N_Vector yscale,   // [in] state scaler
                   N_Vector ydot,      // [in] state residual
                   N_Vector ydotscake,   // [in] state residual scaler
                   void *user_data, // [in/out]
                   N_Vector tmp1, N_Vector tmp2)
{
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorAFSetup(N_Vector y,      // [in]  state vector
                   N_Vector yscale,   // [in] state scaler
                   N_Vector ydot,      // [in] state residual
                   N_Vector ydotscake,   // [in] state residual scaler
                   void *user_data) // [in/out]
{
#endif
  FlameParams *params    = (FlameParams *)user_data;
  const int num_local_points   = params->num_local_points_;
  const int num_states   = params->reactor_->GetNumStates();
  const int num_total_points = params->z_.size();
  const int num_species =  params->reactor_->GetNumSpecies();
  const int num_nonzeros = params->reactor_->GetJacobianSize();
  const int num_states_local = params->num_states_local_;
  double *y_ptr          = NV_DATA_P(y); //_S // caution: assumes realtype == double
  int error_flag = 0;
  MPI_Comm comm = params->comm_;
  int nodeDest;
  int my_pe = params->my_pe_;

  const double constant = params->jacobian_constant_;
  const int nover = params->nover_;

  // Transport Jacobian, evaluated analytically
  for(int j=0; j<num_local_points*5*num_states; j++)
      params->banded_jacobian_[j] = 0.0;

  std::vector<double> dz, inv_dz;
  dz.assign( num_local_points+(2*nover), 0.0);
  inv_dz.assign( num_local_points+(2*nover), 0.0);
  for (int j=0; j<num_local_points+2*nover; ++j) {
    dz[j] = params->dz_local_[j];
    inv_dz[j] = params->inv_dz_local_[j];
  }

  for(int j=0; j<num_local_points; j++) {
    int jext = j + nover;
    int jglobal = j + my_pe*num_local_points;

    double b=0,c=0,d=0; //coefficients of j+1, j, j-1 terms
    double b_m=0,c_m=0,d_m=0;
    double c_p=0,b_p=0,d_p=0;
    // Centered
    b = dz[jext]/dz[jext+1]/(dz[jext]+dz[jext+1]);
    c = (dz[jext+1]-dz[jext])/dz[jext+1]/dz[jext];
    d = -dz[jext+1]/dz[jext]/(dz[jext]+dz[jext+1]);
    b_m = dz[jext-1]/dz[jext]/(dz[jext-1]+dz[jext]);
    c_m = (dz[jext]-dz[jext-1])/dz[jext]/dz[jext-1];
    d_m = -dz[jext]/dz[jext-1]/(dz[jext-1]+dz[jext]);
    b_p = dz[jext+1]/dz[jext+2]/(dz[jext+1]+dz[jext+2]);
    c_p = (dz[jext+2]-dz[jext+1])/dz[jext+2]/dz[jext+1];
    d_p = -dz[jext+2]/dz[jext+1]/(dz[jext+1]+dz[jext+2]);

    // Second derivative
    // Centered three-point stencil for now
    //double bb=0, cc=0, dd=0, bb_m=0, dd_p=0;
    double cc=0, bb_m=0, dd_p=0;
    double denom = dz[jext]*dz[jext+1]*(dz[jext]+dz[jext+1]);
    cc = -2.0*(dz[jext]+dz[jext+1])/denom;
    double denom_m = dz[jext-1]*dz[jext]*(dz[jext-1]+dz[jext]);
    double denom_p = dz[jext+1]*dz[jext+2]*(dz[jext+1]+dz[jext+2]);
    bb_m =  2.0*dz[jext-1]/denom_m;
    dd_p =  2.0*dz[jext+2]/denom_p;

    // Species
    for(int k=0; k<num_species; k++) {
      // Diagonal drhs_j/dY_j
      params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] =
        0.5*cc*params->dissipation_rate_[jext]/
        params->species_lewis_numbers_[jext*num_species + k];

      //Non-unity Le
      if(!params->unity_Lewis_) {
        if(params->full_equations_) {
          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] +=
            0.5*params->dissipation_rate_[jext]/
            params->species_lewis_numbers_[jext*num_species + k]*
            params->mixture_molecular_mass_laplacian_[jext]/
            params->mixture_molecular_mass_[jext];

          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] -=
            0.5*params->dissipation_rate_[jext]/
            params->sum_mass_fraction_over_Lewis_laplacian_[jext];

          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] -=
            0.5*params->dissipation_rate_[jext]/
            params->sum_mass_fraction_over_Lewis_[jext]*
            params->mixture_molecular_mass_laplacian_[jext]/
            params->mixture_molecular_mass_[jext];
        }
        params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] +=
          c*params->convection_velocity_[jext]*
          (1.0/params->species_lewis_numbers_[jext*num_species + k] - 1.0);

        params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] +=
          params->mixture_molecular_mass_grad_[jext]/
          params->species_lewis_numbers_[jext*num_species + k]*(
            params->convection_velocity_[jext]/
            params->mixture_molecular_mass_[jext] +
            0.5*c*params->dissipation_rate_[jext]/
            params->mixture_molecular_mass_[jext]);

        params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] -=
          params->sum_mass_fraction_over_Lewis_grad_[jext]*(
            params->convection_velocity_[jext] +
            0.5*c*params->dissipation_rate_[jext]);

        params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] -=
          params->mixture_molecular_mass_grad_[jext]*(
            params->convection_velocity_[jext]*
            params->sum_mass_fraction_over_Lewis_[jext]/
            params->mixture_molecular_mass_[jext] +
            0.5*c*params->dissipation_rate_[jext]*
            params->sum_mass_fraction_over_Lewis_[jext]/
            params->mixture_molecular_mass_[jext]);

      }

      if(params->soot_) {
        params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] +=
          -params->y_ext_[jext*num_states+num_species]*
          params->rho_dot_[jext] +
          params->y_ext_[jext*num_states+num_species]*
          params->rho_dot_[jext]*params->z_[jglobal]*c;
      }

      // drhs_j-1/dY_j
      if(jglobal > 0) {
        params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 - 1] =
          0.5*bb_m*params->dissipation_rate_[jext-1]/
          params->species_lewis_numbers_[(jext-1)*num_species + k];

        //Non-unity Le
        if(!params->unity_Lewis_) {
          if(params->full_equations_) {
            params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 - 1] +=
              0.5*params->dissipation_rate_[jext-1]/
              params->species_lewis_numbers_[(jext-1)*num_species + k]*
              params->mixture_molecular_mass_laplacian_[jext-1]/
              params->mixture_molecular_mass_[jext-1];

            params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 - 1] -=
              0.5*params->dissipation_rate_[jext-1]/
              params->sum_mass_fraction_over_Lewis_laplacian_[jext-1];

            params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 - 1] -=
              0.5*params->dissipation_rate_[jext-1]/
              params->sum_mass_fraction_over_Lewis_[jext-1]*
              params->mixture_molecular_mass_laplacian_[jext-1]/
              params->mixture_molecular_mass_[jext-1];
          }
          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 - 1] +=
            b_m*params->convection_velocity_[jext-1]*
            (1.0/params->species_lewis_numbers_[(jext-1)*num_species + k] - 1.0);

          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 - 1] +=
            params->mixture_molecular_mass_grad_[jext-1]/
            params->species_lewis_numbers_[(jext-1)*num_species + k]*(
              params->convection_velocity_[jext-1]/
              params->mixture_molecular_mass_[jext-1] +
              0.5*b_m*params->dissipation_rate_[jext-1]/
              params->mixture_molecular_mass_[jext-1]);

          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 - 1] -=
            params->sum_mass_fraction_over_Lewis_grad_[jext-1]*(
              params->convection_velocity_[jext-1] +
              0.5*b_m*params->dissipation_rate_[jext-1]);

          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 - 1] -=
            params->mixture_molecular_mass_grad_[jext-1]*(
              params->convection_velocity_[jext-1]*
              params->sum_mass_fraction_over_Lewis_[jext-1]/
              params->mixture_molecular_mass_[jext-1] +
              0.5*b_m*params->dissipation_rate_[jext-1]*
              params->sum_mass_fraction_over_Lewis_[jext-1]/
              params->mixture_molecular_mass_[jext-1]);

        }

        if(params->soot_) {
          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 - 1] +=
            params->y_ext_[(jext-1)*num_states+num_species]*
            params->rho_dot_[jext-1]*params->z_[jglobal-1]*b_m;
        }

      }

      // drhs_j+1/dY_j
      if(jglobal < num_total_points-1) {
        params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 1] =
          0.5*dd_p*params->dissipation_rate_[jext+1]/
          params->species_lewis_numbers_[(jext+1)*num_species + k];

        //Non-unity Le
        if(!params->unity_Lewis_) {
          if(params->full_equations_) {

            params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 1] +=
              0.5*params->dissipation_rate_[jext+1]/
              params->species_lewis_numbers_[(jext+1)*num_species + k]*
              params->mixture_molecular_mass_laplacian_[jext+1]/
              params->mixture_molecular_mass_[jext+1];

            params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 1] -=
              0.5*params->dissipation_rate_[jext+1]/
              params->sum_mass_fraction_over_Lewis_laplacian_[jext+1];

            params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 1] -=
              0.5*params->dissipation_rate_[jext+1]/
              params->sum_mass_fraction_over_Lewis_[jext+1]*
              params->mixture_molecular_mass_laplacian_[jext+1]/
              params->mixture_molecular_mass_[jext+1];
          }

          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 1] +=
            d_p*params->convection_velocity_[jext+1]*
            (1.0/params->species_lewis_numbers_[(jext+1)*num_species + k] - 1.0);

          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 1] +=
            params->mixture_molecular_mass_grad_[jext+1]/
            params->species_lewis_numbers_[(jext+1)*num_species + k]*(
              params->convection_velocity_[jext+1]/
              params->mixture_molecular_mass_[jext+1] +
              0.5*d_p*params->dissipation_rate_[jext+1]/
              params->mixture_molecular_mass_[jext+1]);

          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 1] -=
            params->sum_mass_fraction_over_Lewis_grad_[jext+1]*(
              params->convection_velocity_[jext+1] +
              0.5*d_p*params->dissipation_rate_[jext+1]);

          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 1] -=
            params->mixture_molecular_mass_grad_[jext+1]*(
              params->convection_velocity_[jext+1]*
              params->sum_mass_fraction_over_Lewis_[jext+1]/
              params->mixture_molecular_mass_[jext+1] +
              0.5*d_p*params->dissipation_rate_[jext+1]*
              params->sum_mass_fraction_over_Lewis_[jext+1]/
              params->mixture_molecular_mass_[jext+1]);

        }

        if(params->soot_) {
          params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 1] +=
            params->y_ext_[(jext+1)*num_states+num_species]*
            params->rho_dot_[jext+1]*params->z_[jglobal+1]*d_p;
        }
      }

    } // for k<num_species

    // Temperature
    if(params->fix_temperature_) {
      params->banded_jacobian_[(num_states-1)*(num_local_points*5) + j*5 + 1 + 2 + 0] = 1.0;
    } else {
      // Diagonal
      params->banded_jacobian_[(num_states-1)*(num_local_points*5) + j*5 + 1 + 2 + 0] =
        0.5*cc*params->dissipation_rate_[jext] +
        (0.5*c*params->dissipation_rate_[jext]/
         params->mixture_specific_heat_[jext]*
         (b*params->mixture_specific_heat_[jext+1] +
          c*params->mixture_specific_heat_[jext] +
          d*params->mixture_specific_heat_[jext-1])) -
        (0.5*c*params->dissipation_rate_[jext]/
         params->mixture_specific_heat_[jext]*
         params->enthalpy_flux_sum_[jext]);

      // drhs_j-1/dY_j
      if(jglobal > 0) {
        params->banded_jacobian_[(num_states-1)*(num_local_points*5) + j*5 + 1 + 2 - 1] =
          0.5*bb_m*params->dissipation_rate_[jext-1] +
          (0.5*b_m*params->dissipation_rate_[jext-1]/
           params->mixture_specific_heat_[jext-1]*
           (b_m*params->mixture_specific_heat_[jext] +
            c_m*params->mixture_specific_heat_[jext-1] +
            d_m*params->mixture_specific_heat_[jext-2])) -
          (0.5*b_m*params->dissipation_rate_[jext-1]/
           params->mixture_specific_heat_[jext-1]*
           params->enthalpy_flux_sum_[jext-1]);
      }

      // drhs_j+1/dY_j
      if(jglobal < num_total_points-1) {
        params->banded_jacobian_[(num_states-1)*(num_local_points*5) + j*5 + 1 + 2 + 1] =
        0.5*dd_p*params->dissipation_rate_[jext+1] +
          (0.5*d_p*params->dissipation_rate_[jext+1]/
           params->mixture_specific_heat_[jext+1]*
           (b_p*params->mixture_specific_heat_[jext+2] +
            c_p*params->mixture_specific_heat_[jext+1] +
            d_p*params->mixture_specific_heat_[jext])) -
          (0.5*d_p*params->dissipation_rate_[jext+1]/
           params->mixture_specific_heat_[jext+1]*
           params->enthalpy_flux_sum_[jext+1]);
      }
    } // if fix_temperature
  } // for j<num_local_points

  // Local chemistry Jacobian (relative volume included here)
  if(params->store_jacobian_) {
    params->saved_jacobian_.assign(num_nonzeros*num_local_points, 0.0);
    // Get Jacobian using Zero-RK
    if(params->fix_temperature_) {
      for(int j=0; j<num_local_points; ++j)
        params->reactor_->GetJacobianDiffusionSteadyFixT(
          &y_ptr[j*num_states],
          &params->step_limiter_[0],
          &params->saved_jacobian_[j*num_nonzeros]);
    } else {
      for(int j=0; j<num_local_points; ++j)
          params->reactor_->GetJacobianDiffusionSteady(
            &y_ptr[j*num_states],
            &params->step_limiter_[0],
            &params->saved_jacobian_[j*num_nonzeros]);
    }
    // Add soot terms
    if(params->soot_) {
      std::vector<double> soot_jacobian;
      soot_jacobian.assign(num_states*num_states, 0.0);

      for(int j=0; j<num_local_points; ++j) {

        ComputeSootJacobian(params,
                            &y_ptr[j*num_states],
                            &soot_jacobian[0]);

        // Add diagonal only for species?
        for(int k=0; k<num_species; ++k) {
          params->saved_jacobian_[j*num_nonzeros+params->diagonal_id_[k]] +=
            soot_jacobian[k*num_states + k];
        }
        // Add relative volume column
        for(int k=0; k<num_states; ++k)
          params->saved_jacobian_[j*num_nonzeros+params->column_sum_[num_species]+k] +=
            soot_jacobian[num_species*num_states + k];

        // Add temperature column
        for(int k=0; k<num_states; ++k)
          params->saved_jacobian_[j*num_nonzeros+params->column_sum_[num_species+1]+k] +=
            soot_jacobian[(num_species+1)*num_states + k];

        // Zero out? shouldn't be necessary
        for(int k=0; k<num_states*num_states; ++k)
          soot_jacobian[k] = 0.0;

      } // for j<num_local_points

    } //if soot

    if(params->pseudo_unsteady_) {
      for(int j=0; j<num_local_points; ++j) {
        // Add -1/dt term to Yi, T
        for(int k=0; k<num_species; ++k) {
          params->saved_jacobian_[j*num_nonzeros+params->diagonal_id_[k]] -= 1.0/params->dt_;
        }
        params->saved_jacobian_[j*num_nonzeros+params->diagonal_id_[num_species+1]] -= 1.0/params->dt_;
      }
    }

    // Add/substract constant on diagonal
    for(int j=0; j<num_local_points; ++j)
      for(int k=0; k<num_states; ++k)
        params->saved_jacobian_[j*num_nonzeros+params->diagonal_id_[k]] -=
          constant;

    // Factorize matrix at each point
    for(int j=0; j<num_local_points; ++j) {
      for(int k=0; k<num_nonzeros; ++k)
        params->reactor_jacobian_[k] = params->saved_jacobian_[j*num_nonzeros+k];

      // factor the numerical jacobian
      if(params->sparse_matrix_[j]->IsFirstFactor()) {
        error_flag =
          params->sparse_matrix_[j]->FactorNewPatternCCS(num_nonzeros,
                                                         &params->row_id_[0],
                                                         &params->column_sum_[0],
                                                         &params->reactor_jacobian_[0]);
      } else {
        error_flag =
          params->sparse_matrix_[j]->FactorSamePattern(
            &params->reactor_jacobian_[0]);
      }
      if(error_flag != 0) {
        params->logger_->PrintF(
          "# DEBUG: At grid point %d (z = %.18g) reactor produced a\n"
          "#        sparse matrix error flag = %d\n",
          j,
          params->z_[j],
          error_flag);

        return error_flag;
      }

    } // for(int j=0; j<num_local_points; ++j)

  } else {
    printf("ERROR: Set store_jacobian to True. Untested otherwise!\n");
    exit(-1);
    // recompute and factor the Jacobian, there is no saved data
    // TODO: offer option for the fake update
    for(int j=0; j<num_local_points; ++j) {
      // Get Jacobian
      params->reactor_->GetJacobianDiffusionSteady(
        &y_ptr[j*num_states],
        &params->step_limiter_[0],
        &params->reactor_jacobian_[0]);

      if(params->pseudo_unsteady_) {
        // Add -1/dt term to Yi, T
        for(int k=0; k<num_species; ++k) {
          params->saved_jacobian_[j*num_nonzeros+params->diagonal_id_[k]] -= 1.0/params->dt_;
        }
        params->saved_jacobian_[j*num_nonzeros+params->diagonal_id_[num_species+1]] -= 1.0/params->dt_;
      }

      //Add/subtract constant on diagonal
      for(int k=0; k<num_states; ++k) {
        params->reactor_jacobian_[params->diagonal_id_[k]] -= constant;
      }

      // factor the numerical jacobian
      if(params->sparse_matrix_[j]->IsFirstFactor()) {
        error_flag =
          params->sparse_matrix_[j]->FactorNewPatternCCS(num_nonzeros,
						 &params->row_id_[0],
						 &params->column_sum_[0],
						 &params->reactor_jacobian_[0]);
      } else {
        error_flag =
          params->sparse_matrix_[j]->FactorSamePattern(
                                                 &params->reactor_jacobian_[0]);
      }
      if(error_flag != 0) {
        params->logger_->PrintF(
          "# DEBUG: At grid point %d (z = %.18g) reactor produced a\n"
          "#        sparse matrix error flag = %d\n",
          j,
          params->z_[j],
          error_flag);

        return error_flag;
      }
    } // for(int j=0; j<num_points; ++j
  } // if(params->store_jacobian_) else
  // Done with chemical Jacobian

  // Add/Subtract identity to/from transport jacobian
  for(int j=0; j<num_local_points; ++j) {
    for(int k=0; k<num_states; ++k) {
      params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] += constant;
    }
  }

  // Multiply by inverse of chemical jacobian
  // TO DO: Make it work when there is no saved jacobian
  double inverse_chem_jacobian;
  for(int j=0; j<num_local_points; ++j) {
    for(int k=0; k<num_states; ++k) {
      inverse_chem_jacobian = 1.0/params->saved_jacobian_[j*num_nonzeros+params->diagonal_id_[k]];
      params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] *= inverse_chem_jacobian;
      params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 - 1] *= inverse_chem_jacobian;
      params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 1] *= inverse_chem_jacobian;
    }
  }

  // Add identity matrix
  for(int j=0; j<num_local_points; ++j) {
    for(int k=0; k<num_states; ++k) {
      params->banded_jacobian_[k*(num_local_points*5) + j*5 + 1 + 2 + 0] += 1.0;
    }
  }

  // Reorganize transport Jacobian for more efficient parallelization prior to factorization
  // Communications to get banded_jacobian2
  long int dsize = num_local_points*5;

  for(int j=0; j<num_states; ++j) {
    nodeDest = j/params->num_states_per_proc_;
    int jlocal = j % params->num_states_per_proc_;
    int start_band = j*(num_local_points*5);
    int start_band2 = jlocal*(num_total_points*5);

    MPI_Gather(&params->banded_jacobian_[start_band],
               dsize,
               PVEC_REAL_MPI_TYPE,
               &params->banded_jacobian2_[start_band2],
               dsize,
               PVEC_REAL_MPI_TYPE,
               nodeDest,
               comm);
  }

  // TODO: Do without jacobian2
  for(int j=0; j<num_states_local; ++j) {
    for(int i=0; i<num_total_points; ++i) {
      for(int s=0; s<4; ++s) {
        params->banded_jacobian_serial_[j*(num_total_points*4) + i*4 + s] =
          params->banded_jacobian2_[j*(num_total_points*5) + i*5 + s + 1];
      }
    }
  }

  // Factorize
  int dim = num_total_points;
  int one = 1;
  int LDAB = 4;
  for(int j=0; j<num_states_local; ++j) {
    dgbtrf_(&dim,
            &dim,
            &one,
            &one,
            &params->banded_jacobian_serial_[j*num_total_points*4],
            &LDAB,
            &params->pivots_serial_[j*num_total_points],
            &error_flag);
  }//for j<num_states_local

  return 0;
}

#if defined SUNDIALS2
int ReactorAFSolve(N_Vector y,      // [in]  state vector
                   N_Vector yscale,   // [in] state scaler
                   N_Vector ydot,      // [in] state residual
                   N_Vector ydotscake,   // [in] state residual scaler
                   N_Vector vv, // [in/out] rhs/solution vecotr
                   void *user_data, // [in/out]
                   N_Vector tmp)
{
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorAFSolve(N_Vector y,      // [in]  state vector
                   N_Vector yscale,   // [in] state scaler
                   N_Vector ydot,      // [in] state residual
                   N_Vector ydotscake,   // [in] state residual scaler
                   N_Vector vv, // [in/out] rhs/solution vecotr
                   void *user_data) // [in/out]
{
#endif
  FlameParams *params = (FlameParams *)user_data;
  const int num_local_points  = params->num_local_points_;
  const int num_states  = params->reactor_->GetNumStates();
  const int num_total_points  = params->z_.size();
  const int num_states_local = params->num_states_local_;
  double *solution    = NV_DATA_P(vv);
  int error_flag = 0;
  int start_id=0;

  // Local sparse chemistry
  for(int j=0; j<num_local_points; ++j) {
    error_flag = params->sparse_matrix_[j]->Solve(&solution[start_id],
                                                  &solution[start_id]);
    start_id += num_states;
    if(error_flag != 0) {
      printf("AFSolve sparse matrix error: %d\n", error_flag);
      return error_flag;
    }
  }

  // Banded transport Jacobian
  // Communications for banded_jacobian2
  MPI_Comm comm = params->comm_;
  long int dsize = num_local_points;
  int nodeDest, nodeFrom;

  std::vector<double> solution_allspecies, solution_species;
  solution_allspecies.assign(num_total_points*num_states_local, 0.0);
  solution_species.assign(num_local_points*num_states, 0.0);

  // Reorder solution vector by species
  for(int j=0; j<num_states; ++j)
    for(int i=0; i<num_local_points; ++i)
      solution_species[j*num_local_points+i] = solution[j+i*num_states];

  // Gather all grid points
  for(int j=0; j<num_states; ++j) {
    nodeDest = j/params->num_states_per_proc_;
    int jlocal = j % params->num_states_per_proc_;
    int start_id = j*num_local_points;
    int start_id2 = jlocal*num_total_points;

    MPI_Gather(&solution_species[start_id],
               dsize,
               PVEC_REAL_MPI_TYPE,
               &solution_allspecies[start_id2],
               dsize,
               PVEC_REAL_MPI_TYPE,
               nodeDest,
               comm);
  }

  // Solve banded matrix
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
      cerr << "Solve banded matrix error: " << error_flag << "\n";
  }

  //Scatter back
  for(int j=0; j<num_states; ++j) {
    nodeFrom = j/params->num_states_per_proc_;
    int jlocal = j % params->num_states_per_proc_;
    int start_id = j*num_local_points;
    int start_id2 = jlocal*num_total_points;

    MPI_Scatter(&solution_allspecies[start_id2],
                dsize,
                PVEC_REAL_MPI_TYPE,
                &solution_species[start_id],
                dsize,
                PVEC_REAL_MPI_TYPE,
                nodeFrom,
                comm);
  }

  for(int j=0; j<num_states; ++j)
    for(int i=0; i<num_local_points; ++i)
      solution[j+i*num_states] = solution_species[j*num_local_points+i];

  return error_flag;
}


// -------------------------------------------------------------------------
// SuperLUDIST block-tridiagonal Jacobian
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
  double *y_ptr          = NV_DATA_P(y);
  double *ydot_ptr       = NV_DATA_P(ydot);
  int error_flag = 0;
  double alpha = 1.0e-6;
  double beta = 1.0e-14;
  double delta;

  int my_pe = params->my_pe_;
  int npes  = params->npes_;
  const int nover=params->nover_;

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
        i1global = max(0, jglobal - jstate - num_states);
        i2global = min(jglobal + (num_states-1 - jstate) + num_states, num_total_states-1);
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

  // Perform parallel communication of jacobian
  MPI_Comm comm = params->comm_;
  MPI_Status status;
  long int dsize_jac_bnd = nover*num_states*width;

  // MPI sendrecv
  int nodeDest = my_pe-1;
  if (nodeDest < 0) nodeDest = npes-1;
  int nodeFrom = my_pe+1;
  if (nodeFrom > npes-1) nodeFrom = 0;
  MPI_Sendrecv(&jac_bnd[nover*num_states*width],
               dsize_jac_bnd, PVEC_REAL_MPI_TYPE, nodeDest, 0,
               &jac_bnd[num_states*(num_local_points+nover)*width],
               dsize_jac_bnd, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);

  nodeDest = my_pe+1;
  if (nodeDest > npes-1) nodeDest = 0;
  nodeFrom = my_pe-1;
  if (nodeFrom < 0) nodeFrom = npes-1;
  MPI_Sendrecv(&jac_bnd[num_states*num_local_points*width],
               dsize_jac_bnd, PVEC_REAL_MPI_TYPE, nodeDest, 0,
               &jac_bnd[0], dsize_jac_bnd, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);

  // Get pattern "manually" for now
  // TODO: find a cleaner way
  int innz=0;
  // here j is the ROW and i is the COLUMN
  for (j=0; j<num_local_states; j++) {
    jglobal = j + my_pe*num_local_states;
    jext = j + nover*num_states;
    int jstate = jglobal % num_states;
    i1global = max(0, jglobal - jstate - num_states);
    i2global = min(jglobal + (num_states-1 - jstate) + num_states, num_total_states-1);
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

  // Factorize
  if(params->sparse_matrix_dist_->IsFirstFactor_dist()) {
    error_flag =
      params->sparse_matrix_dist_->FactorNewPatternCCS_dist(num_nonzeros_loc,
                                                            &params->col_id_[0],
                                                            &params->row_sum_[0],
                                                            &params->reactor_jacobian_dist_[0]);
  } else {
    error_flag =
      params->sparse_matrix_dist_->FactorSamePatternCCS_dist(num_nonzeros_loc,
                                                             &params->col_id_[0],
                                                             &params->row_sum_[0],
                                                             &params->reactor_jacobian_dist_[0]);
  } //if first factor

  return error_flag;
}

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
  double *solution = NV_DATA_P(vv);
  int error_flag = 0;

  error_flag = params->sparse_matrix_dist_->Solve_dist(&solution[0],&solution[0]);

  return error_flag;
}


void ErrorFunction(int error_code,
                   const char *module,
                   const char *function,
                   char *msg,
                   void *user_data)
{
  FlameParams *params = (FlameParams *)user_data;
  MPI_Comm comm = params->comm_;

  printf("Error: %s\n",msg);

  if(error_code != -6) //don't abort if max_iter reached
    MPI_Abort(comm, error_code);
  // Add no abort for step length test?
}

double NonLinearConvectUpwind(double velocity,
                              double y_previous,
                              double y_current,
                              double y_next,
                              double inv_dz_prev,
                              double inv_dz)
{
  if(velocity > 0.0) {
    return (velocity*(y_current - y_previous)*inv_dz_prev);
  } else {
    return (velocity*(y_next - y_current)*inv_dz);
  }
}
