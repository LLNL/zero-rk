#include "compute_lewis.h"

int ComputeLewis(FlameParams &params, double *y_ptr)
{

  const int num_states  = params.reactor_->GetNumStates();
  const int num_species = params.inlet_mass_fractions_.size();

  const double ref_temperature = params.ref_temperature_;

  int transport_error;

  //--------------------------------------------------------------------------
  // Compute the interior heat capacity, conductivity, and species mass fluxes.
  // compute the upstream mid point state for the transport calculations
  for(int k=0; k<num_species; ++k) {

    // mid point mass fractions
    params.transport_input_.mass_fraction_[k] = y_ptr[k];

    // mid point mass fraction gradient
      params.transport_input_.grad_mass_fraction_[k] = 0.0;
  }

  // mid point temperature
  params.transport_input_.temperature_ = ref_temperature*
    y_ptr[num_states-1];

  // mid point temperature gradient
  params.transport_input_.grad_temperature_[0] = 0.0;

  // compute the species mass flux at the upstream mid point
  transport_error = params.trans_->GetSpeciesMassFlux(
    params.transport_input_,
    num_species,
    nullptr,
    nullptr,
    &params.species_mass_flux_[0],
    &params.species_lewis_numbers_[0]);
  if(transport_error != transport::NO_ERROR) {
    return transport_error;
  }

  ofstream myfile;
  myfile.open("Lewis_file");
  for (int k=0; k<num_species; ++k) {
    myfile << params.species_lewis_numbers_[k] << "\n";
  }
  myfile.close();

  return 0;
}
