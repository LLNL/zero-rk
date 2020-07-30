#ifndef FLAME_PARAMS_H_
#define FLAME_PARAMS_H_

#include <string>
#include <vector>

#include <file_utilities.h>

#include <zerork/mechanism.h>
#include <reactor/const_pressure_reactor.h>
#include <transport/mass_transport_factory.h>

#include "UnsteadyFlameIFP.h"

class FlameParams
{
 public:
  explicit FlameParams(const std::string &input_name);
  ~FlameParams();

  UnsteadyFlameIFP *parser_;
  ConstPressureReactor *reactor_;
  transport::MassTransportInterface *trans_;
  zerork::utilities::Logger *logger_;

  zerork::mechanism *mechanism_; // TODO: avoid using a separate mechanism

  std::string input_name_;

  std::vector<int> fuel_species_id_;

  std::vector<int> full_species_id_;

  std::vector<double> inlet_mass_fractions_;

  double pressure_;
  double inlet_temperature_;
  double inlet_molecular_mass_;
  double inlet_relative_volume_;
  double ref_temperature_;

  bool use_equilibrium_;

  // data structures setup by SetMemory()
  transport::MassTransportInput transport_input_;
  std::vector<double> inv_molecular_mass_;     // size = num_species
  std::vector<double> species_specific_heats_; // size = num_species
  std::vector<double> species_mass_flux_;      // size = num_species
  std::vector<double> species_lewis_numbers_;      // size = num_species*num_points
  std::vector<double> thermal_conductivity_;   // size = num_points+1
  std::vector<double> mixture_specific_heat_;  // size = num_points

 private:
  void SetInlet();
  void SetMemory();
};


#endif
