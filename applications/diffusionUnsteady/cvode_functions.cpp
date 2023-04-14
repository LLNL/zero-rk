#include "cvode_functions.h"
#include "flame_params.h"
#include "utilities/math_utilities.h"

extern "C" void dgbtrf_(int* dim1, int* dim2, int* nu, int* nl, double* a, int* lda, int* ipiv, int* info);
extern "C" void dgbtrs_(char *TRANS, int *N, int *NRHS, int* nu, int* nl, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

// upwind convective scheme
static double NonLinearConvectUpwind(double velocity,
                                     double y_previous,
                                     double y_current,
                                     double y_next,
                                     double inv_dz_prev,
                                     double inv_dz);
// Main RHS function
int ConstPressureFlame(realtype t,
		       N_Vector y,
		       N_Vector ydot,
		       void *user_data)
{
  FlameParams *params = (FlameParams *)user_data;
  const int num_local_points = params->num_local_points_;
  const int num_states  = params->reactor_->GetNumStates();
  int Nlocal = num_local_points*num_states;

  // All communications performed in Local
  ConstPressureFlameLocal(Nlocal, t, y, ydot, user_data);

  return 0;
}

// RHS function
int ConstPressureFlameLocal(int nlocal,
			    realtype t,
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
  const int num_species = params->fuel_mass_fractions_.size();
  const int num_local_states = num_local_points*num_states;
  const int nover = params->nover_;
  int my_pe = params->my_pe_;
  int npes  = params->npes_;

  const double ref_temperature = params->ref_temperature_;

  std::vector<double> enthalpies;
  enthalpies.assign(num_species,0.0);

  std::vector<double> rhs;
  rhs.assign(num_local_points*num_states,0.0);

  std::vector<double> mixture_molecular_mass, conductivity_over_cp;
  std::vector<double>  dissipation_rate_times_rho, sum_mass_fraction_over_Lewis;
  std::vector<double> mass_fraction_over_mixture_mass,  mass_fraction_over_mixture_mass_times_sum;
  mixture_molecular_mass.assign(num_local_points+2*nover, 0.0);
  dissipation_rate_times_rho.assign(num_local_points+2*nover, 0.0);
  conductivity_over_cp.assign(num_local_points+2*nover, 0.0);
  sum_mass_fraction_over_Lewis.assign(num_local_points+2*nover, 0.0);
  mass_fraction_over_mixture_mass.assign((num_local_points+2*nover)*num_species, 0.0);
  mass_fraction_over_mixture_mass_times_sum.assign((num_local_points+2*nover)*num_species, 0.0);

  double mass_fraction_sum;
  double relative_volume_j, convection_velocity_j;
  int transport_error;

  double local_max;
  double thermal_diffusivity;

  std::vector<double> rho_dot;
  rho_dot.assign(num_local_points, 0.0);

  // set the derivative to zero
  for(int j=0; j<num_local_states; ++j) {
    ydot_ptr[j] = 0.0;
    rhs[j] = 0.0;
  }

  // compute the constant pressure reactor source term
  // using Zero-RK
  for(int j=0; j<num_local_points; ++j) {
    params->reactor_->GetTimeDerivativeLimiter(t,
                                               &y_ptr[j*num_states],
                                               &params->step_limiter_[0],
                                               &rhs[j*num_states]);
  }

  // compute soot
  if(params->soot_) {
    for(int j=0; j<num_local_points; ++j) {
      CalcRhoDot(params,
                 &y_ptr[j*num_states],
                 rho_dot[j]);

      UpdateProductionRates(params,
                            &y_ptr[j*num_states],
                            &rhs[j*num_states]);
    }
  }

  //--------------------------------------------------------------------------
  // Perform parallel communications
  MPI_Comm comm = params->comm_;
  MPI_Status status;
  long int dsize = num_states*nover;
  std::vector<double> dz, dzm, inv_dz, inv_dzm;

  dz.assign( num_local_points+(2*nover), 0.0);
  dzm.assign( num_local_points+(2*nover), 0.0);
  inv_dz.assign( num_local_points+(2*nover), 0.0);
  inv_dzm.assign( num_local_points+(2*nover), 0.0);

  // Copy y_ptr data into larger arrays
  for (int j=0; j<num_states*num_local_points; ++j) {
    params->y_ext_[num_states*nover + j] = y_ptr[j];
  }

  for (int j=0; j<num_local_points+2*nover; ++j) {
    dz[j] = params->dz_local_[j];
    dzm[j] = params->dzm_local_[j];
    inv_dz[j] = params->inv_dz_local_[j];
    inv_dzm[j] = params->inv_dzm_local_[j];
  }

  // Update ghost cells with send/receive
  int nodeDest = my_pe-1;
  if (nodeDest < 0) nodeDest = npes-1;
  int nodeFrom = my_pe+1;
  if (nodeFrom > npes-1) nodeFrom = 0;
  MPI_Sendrecv(&params->y_ext_[nover*num_states],
               dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
               &params->y_ext_[num_states*(num_local_points+nover)],
               dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);

  nodeDest = my_pe+1;
  if (nodeDest > npes-1) nodeDest = 0;
  nodeFrom = my_pe-1;
  if (nodeFrom < 0) nodeFrom = npes-1;
  MPI_Sendrecv(&params->y_ext_[num_states*num_local_points],
               dsize, PVEC_REAL_MPI_TYPE, nodeDest, 0,
               &params->y_ext_[0],
               dsize, PVEC_REAL_MPI_TYPE, nodeFrom, 0, comm, &status);

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
  for(int j=0; j<num_local_points+2*nover; ++j) {
    int jlocal = j - nover;
    int jext = j;
    int jglobal = jlocal + my_pe*num_local_points;

    relative_volume_j   = params->y_ext_[jext*num_states+num_species];

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
    // only used to get species Lewis numbers in this case
    // TO DO: get Lewis numbers directly
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
    mixture_molecular_mass[j] = 1.0/mass_fraction_weight_sum;

    // compute mixture conductivity over heat capacitiy
    conductivity_over_cp[j] = params->thermal_conductivity_[j]/
      params->mixture_specific_heat_[j];

    // compute dissipation rate multiplied by density
    dissipation_rate_times_rho[j] = params->dissipation_rate_[j]/relative_volume_j;

    // compute sum of mass fraction over Lewis number
    double local_sum_Y_over_Le = 0.0;
    for(int k=0; k<num_species; ++k) {
      local_sum_Y_over_Le += params->y_ext_[jext*num_states+k]/
        params->species_lewis_numbers_[j*num_species + k];
    }
    sum_mass_fraction_over_Lewis[j] = local_sum_Y_over_Le;

    // compute mass fraction over mixture molecular mass and
    // mass fraction over mixture mass times sum of mass fraction over Lewis
    for(int k=0; k<num_species; ++k) {
      mass_fraction_over_mixture_mass[j*num_species+k] = params->y_ext_[jext*num_states+k]/
        mixture_molecular_mass[j];
      mass_fraction_over_mixture_mass_times_sum[j*num_species+k] =
        mass_fraction_over_mixture_mass[j*num_species+k]*sum_mass_fraction_over_Lewis[j];
    }

  } // for j<num_local_points+2*nover

  //--------------------------------------------------------------------------
  // Pre-compute derivatives of sum(Y_i/Le_i) and W
  std::vector<double> sum_mass_fraction_over_Lewis_grad,
    sum_mass_fraction_over_Lewis_laplacian,
    mixture_molecular_mass_grad,
    mixture_molecular_mass_laplacian;

  sum_mass_fraction_over_Lewis_grad.assign(num_local_points, 0.0);
  sum_mass_fraction_over_Lewis_laplacian.assign(num_local_points, 0.0);
  mixture_molecular_mass_grad.assign(num_local_points, 0.0);
  mixture_molecular_mass_laplacian.assign(num_local_points, 0.0);

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

    for(int k=0; k<num_species; ++k){
      sum_mass_fraction_over_Lewis_grad[j] +=
        (b*params->y_ext_[(jext+1)*num_states+k] +
         c*params->y_ext_[jext*num_states+k] +
         d*params->y_ext_[(jext-1)*num_states+k])/
        params->species_lewis_numbers_[j*num_species + k];

      sum_mass_fraction_over_Lewis_laplacian[j] +=
        (bb*params->y_ext_[(jext+1)*num_states+k] +
         cc*params->y_ext_[jext*num_states+k] +
         dd*params->y_ext_[(jext-1)*num_states+k])/
        params->species_lewis_numbers_[j*num_species + k];

    }

    mixture_molecular_mass_grad[j] =
      (b*mixture_molecular_mass[jext+1] +
       c*mixture_molecular_mass[jext] +
       d*mixture_molecular_mass[jext-1]);

    mixture_molecular_mass_laplacian[j] =
      (bb*mixture_molecular_mass[jext+1] +
       cc*mixture_molecular_mass[jext] +
       dd*mixture_molecular_mass[jext-1]);

  }

  //--------------------------------------------------------------------------
  // Compute convective and diffusive terms for species and temperature
  for(int j=0; j<num_local_points; ++j) {
    int jext = j + nover;
    int jlocal = j + nover;
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
      //fix temperature corresponds to modified YSI flamelet equations
      convection_velocity_j = -270.0*pow(params->z_[jglobal],1.8)*
        pow(1.0-params->z_[jglobal],0.75);
    } else {
      convection_velocity_j = 0.25*relative_volume_j*(
        (b*dissipation_rate_times_rho[jlocal+1] +
         c*dissipation_rate_times_rho[jlocal] +
         d*dissipation_rate_times_rho[jlocal-1]) +
        dissipation_rate_times_rho[jlocal]/conductivity_over_cp[jlocal]*
        (b*conductivity_over_cp[jlocal+1] +
         c*conductivity_over_cp[jlocal] +
         d*conductivity_over_cp[jlocal-1]));
    }
    params->convection_velocity_[jext] = convection_velocity_j;

    // compute the species convection/diffusion terms
    for(int k=0; k<num_species; ++k) {

      // Diffusion - 2
      rhs[j*num_states+k] += 0.5*params->dissipation_rate_[jlocal]/
        params->species_lewis_numbers_[jlocal*num_species + k]*
        (bb*params->y_ext_[(jext+1)*num_states+k] +
         cc*params->y_ext_[jext*num_states+k] +
         dd*params->y_ext_[(jext-1)*num_states+k]);

      // Non-unity Lewis terms below
      if(!params->unity_Lewis_) {

        // These three terms can affect convergence.
        // Option to neglect them at your own risk
        if(params->full_equations_) {
        // Diffusion - 4
        rhs[j*num_states+k] += 0.5*params->dissipation_rate_[jlocal]/
          params->species_lewis_numbers_[jlocal*num_species + k]*
          params->y_ext_[jext*num_states+k]/
          mixture_molecular_mass[jlocal]*
          mixture_molecular_mass_laplacian[j];

        // Diffusion correction - 5.1
        rhs[j*num_states+k] -= 0.5*params->dissipation_rate_[jlocal]*
          params->y_ext_[jext*num_states+k]*
          sum_mass_fraction_over_Lewis_laplacian[j];

        // Diffusion correction - 5.2
        rhs[j*num_states+k] -= 0.5*params->dissipation_rate_[jlocal]*
          params->y_ext_[jext*num_states+k]/
          mixture_molecular_mass[jlocal]*
          sum_mass_fraction_over_Lewis[jlocal]*
          mixture_molecular_mass_laplacian[j];
        }

        // Mass fraction convection velocity - 6
        double uConv = convection_velocity_j*
          (1.0/params->species_lewis_numbers_[jlocal*num_species + k] - 1.0);

        rhs[j*num_states+k] -= NonLinearConvectUpwind(-uConv,
                                                      params->y_ext_[(jext-1)*num_states+k],
                                                      params->y_ext_[jext*num_states+k],
                                                      params->y_ext_[(jext+1)*num_states+k],
                                                      inv_dz[jext],
                                                      inv_dz[jext+1]);

        // Molar mass convection velocity - 7
        rhs[j*num_states+k] += mixture_molecular_mass_grad[j]/
           params->species_lewis_numbers_[jlocal*num_species + k]*(
             convection_velocity_j*params->y_ext_[jext*num_states+k]/
             mixture_molecular_mass[jlocal] +
             0.5*params->dissipation_rate_[jlocal]*
             (b*mass_fraction_over_mixture_mass[(jlocal+1)*num_species+k] +
              c*mass_fraction_over_mixture_mass[(jlocal)*num_species+k] +
              d*mass_fraction_over_mixture_mass[(jlocal-1)*num_species+k]));

        // Mass fraction convection velocity - 8
        rhs[j*num_states+k] -= sum_mass_fraction_over_Lewis_grad[j]*(
            convection_velocity_j*params->y_ext_[jext*num_states+k] +
            0.5*params->dissipation_rate_[jlocal]*
            (a*params->y_ext_[(jext+2)*num_states+k] +
             b*params->y_ext_[(jext+1)*num_states+k] +
             c*params->y_ext_[jext*num_states+k] +
             d*params->y_ext_[(jext-1)*num_states+k] +
             e*params->y_ext_[(jext-2)*num_states+k]));

        // Molar mass convection correction - 9
        rhs[j*num_states+k] -= mixture_molecular_mass_grad[j]*(
             convection_velocity_j*params->y_ext_[jext*num_states+k]/
             mixture_molecular_mass[jlocal]*sum_mass_fraction_over_Lewis[jlocal] +
             0.5*params->dissipation_rate_[jlocal]*
             (b*mass_fraction_over_mixture_mass_times_sum[(jlocal+1)*num_species+k] +
              c*mass_fraction_over_mixture_mass_times_sum[(jlocal)*num_species+k] +
              d*mass_fraction_over_mixture_mass_times_sum[(jlocal-1)*num_species+k]));

      } // if !unity_Lewis

      if(params->soot_) {
        // -rhodot*(Y_k - Z*dY_k/dZ)
        rhs[j*num_states+k] -= rho_dot[j]*relative_volume_j*
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
      rhs[(j+1)*num_states-1] = 0.0;
    } else {
      // Heat conduction - 2
      rhs[(j+1)*num_states-1] += 0.5*params->dissipation_rate_[jlocal]*
        (bb*params->y_ext_[(jext+2)*num_states-1] +
         cc*params->y_ext_[(jext+1)*num_states-1] +
         dd*params->y_ext_[jext*num_states-1]);

      // Heat condution - 3
      rhs[(j+1)*num_states-1] += 0.5*params->dissipation_rate_[jlocal]/
        params->mixture_specific_heat_[jlocal]* //divide by Cp
        (b*params->mixture_specific_heat_[jlocal+1] +
         c*params->mixture_specific_heat_[jlocal] +
         d*params->mixture_specific_heat_[jlocal-1])*
        (b*params->y_ext_[(jext+2)*num_states-1] +
         c*params->y_ext_[(jext+1)*num_states-1] +
         d*params->y_ext_[jext*num_states-1]);

      // Enthalpy flux - 7
      double enthalpy_flux_sum = 0.0;
      for(int k=0; k<num_species; k++) {

        enthalpy_flux_sum -= params->species_specific_heats_[num_species*jlocal + k]/
          params->species_lewis_numbers_[jlocal*num_species + k]*
          (a*params->y_ext_[(jext+2)*num_states+k] +
           b*params->y_ext_[(jext+1)*num_states+k] +
           c*params->y_ext_[jext*num_states+k] +
           d*params->y_ext_[(jext-1)*num_states+k] +
           e*params->y_ext_[(jext-2)*num_states+k]);

        if(!params->unity_Lewis_) {

          enthalpy_flux_sum += params->mixture_specific_heat_[jlocal]/
            params->species_lewis_numbers_[jlocal*num_species + k]*
            (a*params->y_ext_[(jext+2)*num_states+k] +
             b*params->y_ext_[(jext+1)*num_states+k] +
             c*params->y_ext_[jext*num_states+k] +
             d*params->y_ext_[(jext-1)*num_states+k] +
             e*params->y_ext_[(jext-2)*num_states+k]);

          enthalpy_flux_sum += (params->mixture_specific_heat_[jlocal] -
                                params->species_specific_heats_[jlocal*num_species + k])/
            params->species_lewis_numbers_[jlocal*num_species + k]*
            params->y_ext_[jext*num_states+k]/
            mixture_molecular_mass[jlocal]*
            (b*mixture_molecular_mass[jlocal+1] +
             c*mixture_molecular_mass[jlocal] +
             d*mixture_molecular_mass[jlocal-1]);

        } // if !unity_Lewis
      } // for k<num_species

      // Save enthalpy_flux_sum for use in jacobian
      params->enthalpy_flux_sum_[jext] = enthalpy_flux_sum;

      rhs[(j+1)*num_states-1] -= 0.5*params->dissipation_rate_[jlocal]/
        params->mixture_specific_heat_[jlocal]*
        enthalpy_flux_sum*(b*params->y_ext_[(jext+2)*num_states-1] +
                           c*params->y_ext_[(jext+1)*num_states-1] +
                           d*params->y_ext_[jext*num_states-1]);

      if(params->soot_) {
        // rhodot*Z*dcpT/dZ (divided by rho*cp for correct units)
        rhs[(j+1)*num_states-1] += rho_dot[j]*params->z_[jglobal]*relative_volume_j/
          params->mixture_specific_heat_[jlocal]*
          (b*params->y_ext_[(jext+2)*num_states-1]*params->mixture_specific_heat_[jlocal+1] +
           c*params->y_ext_[(jext+1)*num_states-1]*params->mixture_specific_heat_[jlocal] +
           d*params->y_ext_[jext*num_states-1]*params->mixture_specific_heat_[jlocal-1]);
      }

    } // if else fix_temperature

  } // for(int j=1; j<num_local_points; ++j) // loop computing rhs

  // -------------------------------------------------------------------------
  // Compute the rate of change of the relative volume using the ideal
  // gas equation of state, and the FINAL derivatives of temperature and
  // mass fractions.
  //
  // dv/dt = v/T * dT/dt + RuT/p * \sum_i (1/mw[i] * dy[i]/dt)
  // dv/dt   current units [m^3/kg/s]
  const double RuTref_p = params->reactor_->GetGasConstant()*
    params->ref_temperature_/params->parser_->pressure();

  for(int j=0; j<num_local_points; ++j) {
    int rvol_id = j*num_states+num_species; // relative volume index of pt j
    int temp_id = rvol_id+1;                // temperature index of pt j

    mass_fraction_sum = 0.0;
    for(int k=0; k<num_species; ++k) {
      ydot_ptr[j*num_states+k] = rhs[j*num_states+k];
      mass_fraction_sum += params->inv_molecular_mass_[k]*ydot_ptr[j*num_states+k];
    }
    ydot_ptr[temp_id] = rhs[temp_id];

    ydot_ptr[rvol_id] = y_ptr[rvol_id]*ydot_ptr[temp_id]/y_ptr[temp_id] +
      RuTref_p*y_ptr[temp_id]*mass_fraction_sum;
  }

  // -------------------------------------------------------------------------
  // Parallel communication for finite difference Jacobian
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

  // -------------------------------------------------------------------------
  // Parallel communication of enthalpy flux sum for analytical transport Jacobian
  if(params->integrator_type_ == 3 && params->implicit_transport_) {
    MPI_Status status;
    long int dsize = nover;

    if(npes>1) {
      // MPI sendrecv
      int nodeDest = my_pe-1;
      if (nodeDest < 0) nodeDest = npes-1;
      int nodeFrom = my_pe+1;
      if (nodeFrom > npes-1) nodeFrom = 0;
      MPI_Sendrecv(&params->enthalpy_flux_sum_[nover],
                   dsize, PVEC_REAL_MPI_TYPE,
                   nodeDest, 0,
                   &params->enthalpy_flux_sum_[(num_local_points+nover)],
                   dsize, PVEC_REAL_MPI_TYPE,
                   nodeFrom, 0, comm, &status);

      nodeDest = my_pe+1;
      if (nodeDest > npes-1) nodeDest = 0;
      nodeFrom = my_pe-1;
      if (nodeFrom < 0) nodeFrom = npes-1;
      MPI_Sendrecv(&params->enthalpy_flux_sum_[num_local_points],
                   dsize, PVEC_REAL_MPI_TYPE,
                   nodeDest, 0,
                   &params->enthalpy_flux_sum_[0],
                   dsize, PVEC_REAL_MPI_TYPE,
                   nodeFrom, 0, comm, &status);
    }
  }

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

#if defined SUNDIALS2
int ReactorPreconditionerSetup(realtype t,      // [in] ODE system time
                               N_Vector y,      // [in] ODE state vector
                               N_Vector ydot,   // [in] ODE state derivative
                               booleantype jok,
			       booleantype *new_j,
			       realtype gamma,
		               void *user_data,
                               N_Vector tmp1,
                               N_Vector tmp2,
                               N_Vector tmp3) // [in/out]
{
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorPreconditionerSetup(realtype t,      // [in] ODE system time
                               N_Vector y,      // [in] ODE state vector
                               N_Vector ydot,   // [in] ODE state derivative
                               booleantype jok,
			       booleantype *new_j,
			       realtype gamma,
		               void *user_data) // [in/out]
{
#endif
  FlameParams *params    = (FlameParams *)user_data;
  const int num_local_points   = params->num_local_points_;
  const int num_states   = params->reactor_->GetNumStates();
  const int num_total_points = params->z_.size();
  const int num_species = params->fuel_mass_fractions_.size();
  const int num_nonzeros = params->reactor_->GetJacobianSize();
  const int num_states_local = params->num_states_local_;
  double *y_ptr          = NV_DATA_P(y); //_S // caution: assumes realtype == double
  int error_flag = 0;
  int my_pe = params->my_pe_;

  if(params->store_jacobian_) {

    if(!jok) {
      // The Jacobian is not okay, need to recompute
      params->saved_jacobian_.assign(num_nonzeros*num_local_points, 0.0);
      for(int j=0; j<num_local_points; ++j) {
        params->reactor_->GetJacobianLimiter(t,
					     &y_ptr[j*num_states],
					     &params->step_limiter_[0],
					     &params->saved_jacobian_[j*num_nonzeros]);
      }
      (*new_j) = true;
    } else {
      (*new_j) = false;
    }

    for(int j=0; j<num_local_points; ++j) {

      // compute I - gamma*J
      // this could be updated with blas routines
      for(int k=0; k<num_nonzeros; ++k) {
        params->reactor_jacobian_[k] =
          -gamma*params->saved_jacobian_[j*num_nonzeros+k];
      }
      for(int k=0; k<num_states; ++k) {
        params->reactor_jacobian_[params->diagonal_id_[k]] += 1.0;
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
          "# DEBUG: At t = %.18g [s],\n"
          "#        grid point %d (z = %.18g [m]) reactor produced a\n"
          "#        sparse matrix error flag = %d\n",
          t,
          j,
          params->z_[j],
          error_flag);

        return error_flag;
      }

     } // for(int j=0; j<num_local_points; ++j)

  } else {
    // recompute and factor the Jacobian, there is no saved data
    // TODO: offer option for the fake update
    for(int j=0; j<num_local_points; ++j) {
      params->reactor_->GetJacobianLimiter(t,
					   &y_ptr[j*num_states],
					   &params->step_limiter_[0],
					   &params->reactor_jacobian_[0]);

      // compute I - gamma*J
      // this could be updated with blas routines
      for(int k=0; k<num_nonzeros; ++k) {
        params->reactor_jacobian_[k] *= -gamma;
      }
      for(int k=0; k<num_states; ++k) {
        params->reactor_jacobian_[params->diagonal_id_[k]] += 1.0;
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
          "# DEBUG: At t = %.18g [s],\n"
          "#        grid point %d (z = %.18g [m]) reactor produced a\n"
          "#        sparse matrix error flag = %d\n",
          t,
          j,
          params->z_[j],
          error_flag);

        return error_flag;
      }

     } // for(int j=0; j<num_points; ++j)

    (*new_j) = true; // without saving, it is always a new Jacobian

  } // if(params->store_jacobian_) else

  if(params->implicit_transport_) {
    // Fill the jacobian matrix
    int mkeep = params->num_off_diagonals_;
    int storage = 4*mkeep + 1;
    int fillin = 2*mkeep;
    for(int j=0; j<num_local_points*storage*num_states; j++)
      params->banded_jacobian_[j] = 0.0;

    const int nover = params->nover_;
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

    for(int j=0; j<num_local_points; j++) {
      int jlocal = j + nover;
      int jext = j + nover;
      int jglobal = j + my_pe*num_local_points;

      double b=0,c=0,d=0; //coefficients of j+1, j, j-1 terms
      // Centered
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

      // Species
      for(int k=0; k<num_species; k++) {
	// Diagonal drhs_j/dY_j
	params->banded_jacobian_[k*(num_local_points*storage) + j*storage + mkeep + fillin + 0] =
	  1.0 - gamma*(0.5*cc*params->dissipation_rate_[jlocal]/
                       params->species_lewis_numbers_[jlocal*num_species + k]);//add non-unity Le?

	// drhs_j/dY_j-1
	if(jglobal > 0) {
          params->banded_jacobian_[k*(num_local_points*storage) + j*storage + mkeep + fillin + 1] =
	    -gamma*(0.5*dd*params->dissipation_rate_[jlocal]/
                    params->species_lewis_numbers_[(jlocal)*num_species + k]);//add non-unity Le?
	}

	// drhs_j/dY_j+1
	if(jglobal < num_total_points-1) {
          params->banded_jacobian_[k*(num_local_points*storage) + j*storage + mkeep + fillin - 1] =
	    -gamma*(0.5*bb*params->dissipation_rate_[jlocal]/
                    params->species_lewis_numbers_[(jlocal)*num_species + k]);//add non-unity Le?
	}

      } // for k<num_species

      // Temperature
      // Diagonal
      params->banded_jacobian_[(num_species+1)*(num_local_points*storage) + j*storage + mkeep + fillin + 0] =
	1.0 - gamma*(0.5*cc*params->dissipation_rate_[jlocal]+
                     (0.5*c*params->dissipation_rate_[jlocal]/
                      params->mixture_specific_heat_[jlocal]*
                      (b*params->mixture_specific_heat_[jlocal+1] +
                       c*params->mixture_specific_heat_[jlocal] +
                       d*params->mixture_specific_heat_[jlocal-1])) -
                     (0.5*c*params->dissipation_rate_[jlocal]/
                      params->mixture_specific_heat_[jlocal]*
                      params->enthalpy_flux_sum_[jext]) );

      // drhs_j/dY_j-1
      if(jglobal > 0) {
        params->banded_jacobian_[(num_species+1)*(num_local_points*storage) + j*storage + mkeep + fillin + 1] =
	  -gamma*(0.5*dd*params->dissipation_rate_[jlocal] +
                  (0.5*d*params->dissipation_rate_[jlocal]/
                   params->mixture_specific_heat_[jlocal]*
                   (b*params->mixture_specific_heat_[jlocal+1] +
                    c*params->mixture_specific_heat_[jlocal] +
                    d*params->mixture_specific_heat_[jlocal-1])) -
                  (0.5*d*params->dissipation_rate_[jlocal]/
                   params->mixture_specific_heat_[jlocal]*
                   params->enthalpy_flux_sum_[jext]) );
      }

      // drhs_j/dY_j+1
      if(jglobal < num_total_points-1) {
        params->banded_jacobian_[(num_species+1)*(num_local_points*storage) + j*storage + mkeep + fillin - 1] =
	  -gamma*(0.5*bb*params->dissipation_rate_[jlocal] +
                  (0.5*b*params->dissipation_rate_[jlocal]/
                   params->mixture_specific_heat_[jlocal]*
                   (b*params->mixture_specific_heat_[jlocal+1] +
                    c*params->mixture_specific_heat_[jlocal] +
                    d*params->mixture_specific_heat_[jlocal-1])) -
                  (0.5*b*params->dissipation_rate_[jlocal]/
                   params->mixture_specific_heat_[jlocal]*
                   params->enthalpy_flux_sum_[jext]) );
      }

      // Relative volume
      // Diagonal
      params->banded_jacobian_[(num_species)*(num_local_points*storage) + j*storage + mkeep + fillin + 0] = 1.0;
    } // for j<num_local_points

    // Communications to get banded_jacobian2
    MPI_Comm comm = params->comm_;
    long int dsize = num_local_points*5;
    int nodeDest;

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


  } // if(params->implicit_transport_)

  return 0;
}

#if defined SUNDIALS2
int ReactorPreconditionerSolve(realtype t,      // [in] ODE system time
                               N_Vector y,      // [in] ODE state vector
                               N_Vector ydot,   // [in] ODE state derivative
                               N_Vector r,      // [in] jacobian rhs
                               N_Vector z,      // [out]
			       realtype gamma,
			       realtype delta,
			       int lr,
		               void *user_data,
                               N_Vector tmp)    // [in/out]
{
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorPreconditionerSolve(realtype t,      // [in] ODE system time
                               N_Vector y,      // [in] ODE state vector
                               N_Vector ydot,   // [in] ODE state derivative
                               N_Vector r,      // [in] jacobian rhs
                               N_Vector z,      // [out]
			       realtype gamma,
			       realtype delta,
			       int lr,
		               void *user_data)    // [in/out]
{
#endif
  FlameParams *params = (FlameParams *)user_data;
  //  const int num_points  = params->z_.size();
  const int num_local_points  = params->num_local_points_;
  const int num_states  = params->reactor_->GetNumStates();
  const int num_total_points  = params->z_.size();
  const int num_states_local = params->num_states_local_;
  double *rhs         = NV_DATA_P(r);  // pointers to data array for N_Vector
  double *solution    = NV_DATA_P(z);  // pointers to data array for N_Vector
  int error_flag = 0;
  int start_id=0;

  for(int j=0; j<num_local_points; ++j) {
    error_flag = params->sparse_matrix_[j]->Solve(&rhs[start_id],
                                                  &solution[start_id]);
    start_id += num_states;
    if(error_flag != 0) {
      return error_flag;
    }
  }

  if(params->implicit_transport_) {

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

  } // if(params->implicit_transport_)
  return error_flag;
}

// -------------------------------------------------------------------------
// SuperLUDIST block-tridiagonal Jacobian
#if defined SUNDIALS2
int ReactorBBDSetup(realtype t,      // [in] ODE system time
                    N_Vector y,      // [in] ODE state vector
                    N_Vector ydot,   // [in] ODE state derivative
                    booleantype jok,
                    booleantype *new_j,
                    realtype gamma,
                    void *user_data,
                    N_Vector tmp1,
                    N_Vector tmp2,
                    N_Vector tmp3) // [in/out]
{
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorBBDSetup(realtype t,      // [in] ODE system time
                    N_Vector y,      // [in] ODE state vector
                    N_Vector ydot,   // [in] ODE state derivative
                    booleantype jok,
                    booleantype *new_j,
                    realtype gamma,
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
  ConstPressureFlame(t, y, ydot, user_data);

  // Save copy of state vector and rhs
  for (int j=0; j<num_local_states; ++j)
    y_saved[j] = y_ptr[j];

  for (int j=0; j<num_states*(num_local_points+2*nover); ++j)
    rhs_ext_saved[j] = params->rhs_ext_[j];

  // Banded Jacobian
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
    ConstPressureFlame(t, y, ydot, user_data);

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
          params->reactor_jacobian_dist_[innz] = -gamma*jac_bnd[iext*width + jglobal-i+mkeep];
          if(jstate==istate) {
            params->reactor_jacobian_dist_[innz] += 1.0;
          }
          innz++;
        }

      }
      //Off-diagonal blocks
      if (i<jglobal-jstate || i>jglobal+num_states-1-jstate) {
        if (params->dense_to_sparse_offdiag_[dense_id] == 1) {
          params->reactor_jacobian_dist_[innz] = -gamma*jac_bnd[iext*width + jglobal-i+mkeep];
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
int ReactorBBDSolve(realtype t,      // [in] ODE system time
                    N_Vector y,      // [in] ODE state vector
                    N_Vector ydot,   // [in] ODE state derivative
                    N_Vector r,      // [in] jacobian rhs
                    N_Vector z,      // [out]
                    realtype gamma,
                    realtype delta,
                    int lr,
                    void *user_data,
                    N_Vector tmp)    // [in/out]
{
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorBBDSolve(realtype t,      // [in] ODE system time
                    N_Vector y,      // [in] ODE state vector
                    N_Vector ydot,   // [in] ODE state derivative
                    N_Vector r,      // [in] jacobian rhs
                    N_Vector z,      // [out]
                    realtype gamma,
                    realtype delta,
                    int lr,
                    void *user_data)    // [in/out]
{
#endif
  FlameParams *params = (FlameParams *)user_data;
  double *rhs = NV_DATA_P(r);
  double *solution = NV_DATA_P(z);
  int error_flag = 0;

  error_flag = params->sparse_matrix_dist_->Solve_dist(&rhs[0],&solution[0]);

  return error_flag;

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
