#ifndef CONSTANT_LEWIS_H_
#define CONSTANT_LEWIS_H_

#include <zerork/mechanism.h>

#include "mass_transport_factory.h" // abstract base class

namespace transport
{

class ConstantLewis : public MassTransportInterface
{
 public:
  ConstantLewis();
  ~ConstantLewis();
  int Initialize(const std::vector<std::string> &input_files,
                 const std::string &log_name,
                 const double conductivity_multiplier = 1.0);

  int Initialize(zerork::mechanism* mechanism,
                 const std::vector<std::string> &input_files,
                 const std::string &log_name,
                 const double conductivity_multiplier = 1.0);

  int GetMixtureViscosity(const MassTransportInput &input,
                          double *viscosity) const;
  int GetSpeciesViscosity(const MassTransportInput &input,
                          double *viscosity) const;

  int GetMixtureConductivity(const MassTransportInput &input,
                             double *conductivity) const;

  int GetSpeciesConductivity(const MassTransportInput &input,
                             double *conductivity) const;

  int GetSpeciesMassFlux(const MassTransportInput &input,
                         const size_t ld_species_mass_flux,
			 double *conductivity_mix,
			 double *specific_heat_mix,
                         double *species_mass_flux,
			 double *species_lewis_numbers ) const;

  int GetSpeciesMassFluxFrozenThermo(const MassTransportInput &input,
				     const size_t ld_species_mass_flux,
				     double *conductivity_mix,
				     double *specific_heat_mix,
				     double *species_mass_flux,
				     double *species_lewis_numbers ) const;

 private:
  bool initialized_;
  int num_species_;
  std::string log_name_;
  std::vector<double> molecular_mass_;
  std::vector<double> inv_molecular_mass_;
  std::vector<double> lewis_number_;

  std::vector<int> shape_;
  std::vector<double> kOverEps_;
  std::vector<double> sigma_;
  std::vector<double> mu_;
  std::vector<double> alpha_;
  std::vector<double> Zrot_;

  std::vector<double> mucoeff_;

  mutable std::vector<double> species_workspace_;

  double multiplier_;

  zerork::mechanism *mechanism_; // TODO: avoid using a separate mechanism
                                 //       instantiation
  bool mechanism_owner_;
  int ParseTransportFile(const std::string &transport_file,
                         const std::string &ignore_chars,
                         std::string *error_message);

  double omega_D (double t) const;
  double omega_mu (double t) const;
};


} // namespace transport

#endif
