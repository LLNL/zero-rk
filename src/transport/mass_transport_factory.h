#ifndef MASS_TRANSPORT_FACTORY_H_
#define MASS_TRANSPORT_FACTORY_H_

#include <vector>
#include <string>

#include "zerork/mechanism.h"

// Abstract base class conforming to the Google style rules for a pure
// interface:
//
//   1. It has only public pure virtual ("= 0") methods and static methods
//      (but see below for destructor).
//   2. It may not have non-static data members.
//   3. It need not have any constructors defined. If a constructor is
//      provided, it must take no arguments and it must be protected.
//   4. If it is a subclass, it may only be derived from classes that satisfy
//      these conditions and are tagged with the Interface suffix.
//
// An interface class can never be directly instantiated because of the pure
// virtual method(s) it declares. To make sure all implementations of the
// interface can be destroyed correctly, the interface must also declare a
// virtual destructor (in an exception to the first rule, this should not be
// pure). See Stroustrup, The C++ Programming Language, 3rd edition,
// section 12.4 for details.
//
// Note that only clases that satisfy the conditions for a pure interface
// can end with an Interface suffix.

namespace transport {

// Note that not every input parameter will always need to be set for each
// interface
typedef struct
{
  size_t num_dimensions_;
  size_t ld_grad_temperature_;
  size_t ld_grad_pressure_;
  size_t ld_grad_mass_fraction_;
  double temperature_;
  double pressure_;
  double *mass_fraction_;
  double *grad_temperature_;
  double *grad_pressure_;
  double *grad_mass_fraction_;
} MassTransportInput;

// define mass transport interface error codes
const int NO_ERROR         =  0;
const int UNINITIALIZED    = -1;
const int LOG_FILE_ERROR   = -2;
const int INPUT_FILE_ERROR = -3;
const int NON_POSITIVE_TEMPERATURE = -4;

class MassTransportInterface
{
 public:
  virtual ~MassTransportInterface() {};

  virtual int Initialize(const std::vector<std::string> &input_files,
                         const std::string &log_name,
                         const double conductivity_multiplier = 1.0) = 0;

  virtual int Initialize(zerork::mechanism* mechanism,
                         const std::vector<std::string> &input_files,
                         const std::string &log_name,
                         const double conductivity_multiplier = 1.0) = 0;

  virtual int GetMixtureViscosity(const MassTransportInput &input,
                                  double *viscosity) const = 0;
  virtual int GetSpeciesViscosity(const MassTransportInput &input,
                                  double *viscosity) const = 0;

  virtual int GetMixtureConductivity(const MassTransportInput &input,
                                     double *conductivity) const = 0;
  virtual int GetSpeciesConductivity(const MassTransportInput &input,
                                     double *conductivity) const = 0;

  virtual int GetSpeciesMassFlux(const MassTransportInput &input,
                                 const size_t ld_species_mass_flux,
				 double *conductivity_mix,
				 double *specific_heat_mix,
                                 double *species_mass_flux,
				 double *species_lewis_numbers) const = 0;

  virtual int GetSpeciesMassFluxFrozenThermo(const MassTransportInput &input,
					     const size_t ld_species_mass_flux,
					     double *conductivity_mix,
					     double *specific_heat_mix,
					     double *species_mass_flux,
					     double *species_lewis_numbers) const = 0;

};


class InterfaceFactory
{
 public:
  static MassTransportInterface *CreateMassBased(const std::string &type);
};


} // end namespace transport

#endif
