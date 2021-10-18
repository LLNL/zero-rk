#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
//#include <utilities/string_utilities.h>
//#include <utilities/math_utilities.h>

#include "set_initial_conditions.h"

// Set all grid points to the inlet conditions including temperature
void SetConstantInlet(FlameParams &flame_params, double *y)
{
  const int num_species = flame_params.inlet_mass_fractions_.size();

  const double pressure = flame_params.parser_->pressure();
  const double inlet_temperature = flame_params.parser_->inlet_temperature();
  const double ref_temperature = flame_params.parser_->ref_temperature();

  double relative_volume;
  double mixture_molecular_mass = 0.0;
  std::vector<double> species_molecular_mass;

  species_molecular_mass.assign(num_species, 0.0);
  flame_params.reactor_->GetSpeciesMolecularWeight(&species_molecular_mass[0]);
  for(int k=0; k<num_species; ++k) {
    mixture_molecular_mass +=
      flame_params.inlet_mass_fractions_[k]/species_molecular_mass[k];
  }
  mixture_molecular_mass = 1.0/mixture_molecular_mass;
  relative_volume = flame_params.reactor_->GetGasConstant()*inlet_temperature/
    (pressure*mixture_molecular_mass);

  for(int k=0; k<num_species; ++k) {
    y[k] = flame_params.inlet_mass_fractions_[k];
  }
  y[num_species]   = relative_volume;
  y[num_species+1] = inlet_temperature/ref_temperature;

}
