#ifndef FLEXIBLE_TRANSPORT_H_
#define FLEXIBLE_TRANSPORT_H_

#include <zerork/mechanism.h>

#include "mass_transport_factory.h" // abstract base class

namespace transport
{

class FlexibleTransport : public MassTransportInterface
{
 public:
  FlexibleTransport();
  ~FlexibleTransport();
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
                         double *species_lewis_numbers) const;

  int GetSpeciesMassFluxFrozenThermo(const MassTransportInput &input,
                                     const size_t ld_species_mass_flux,
                                     double *conductivity_mix,
                                     double *specific_heat_mix,
                                     double *species_mass_flux,
                                     double *species_lewis_numbers) const;

  void SetMixAvg(bool setting) {mix_avg_ = setting;};
  void SetSoret(bool setting) {soret_ = setting;};

 private:
  int GetSpeciesMassFluxInternal(const MassTransportInput &input,
                                 const size_t ld_species_mass_flux,
                                 double *conductivity_mix,
                                 double *specific_heat_mix,
                                 double *species_mass_flux,
                                 double *species_lewis_numbers,
                                 bool frozen) const;

  bool initialized_;
  bool mix_avg_;
  bool soret_;
  int num_species_;
  std::string log_name_;
  std::vector<double> molecular_mass_;
  std::vector<double> inv_molecular_mass_;

  std::vector<int> shape_;
  std::vector<double> kOverEps_;
  std::vector<double> sqrtkOverEps_;
  std::vector<double> sigma_;
  std::vector<double> mu_;
  std::vector<double> alpha_;
  std::vector<double> Zrot_;

  bool precompute_matrix_terms_;
  std::vector<double> sqrtmass_, diam2_;
  std::vector<double> sqrt2mass_, inv_sqrt1mass_, inv_sum_mass_;

  mutable std::vector<double> lewis_numbers_, Dmass, invDij;
  mutable std::vector<double> DTherm, DeltaI, sqrtmu, inv_sqrtmu;

  std::vector<double> mucoeff_;

  double multiplier_;

  mutable std::vector<double> species_workspace_;

  zerork::mechanism *mechanism_;
  bool mechanism_owner_;
  int ParseTransportFile(const std::string &transport_file,
                         const std::string &ignore_chars,
                         std::string *error_message);

  double omega_D (double t) const;
  double omega_C (double t) const;
  double omega_mu (double t) const;

};

} // namespace transport

#endif
