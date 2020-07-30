#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
#include <utilities/string_utilities.h>
#include <utilities/math_utilities.h>

#include "flame_params.h"


static double NormalizeComposition(const size_t num_elements,
                                   double composition[]);

FlameParams::FlameParams(const std::string &input_name)
{

  int error_code;

  parser_  = NULL;
  reactor_ = NULL;
  trans_   = NULL;
  logger_  = NULL;

  if(!zerork::utilities::FileIsReadable(input_name)) {
    printf("# ERROR: Input file %s is not readable\n",input_name.c_str());
    exit(-1);
  }

  // setup parser
  input_name_ = input_name;
  parser_ = new UnsteadyFlameIFP(input_name);
  if(parser_ == NULL) {
    printf("# ERROR: Parser for file %s is not created\n",input_name.c_str());
    exit(-1);
  }


  // setup constant pressure reactor
  reactor_ = new ConstPressureReactor(parser_->mech_file().c_str(),
                                      parser_->therm_file().c_str(),
                                      parser_->log_file().c_str(),
                                      COMPRESSED_COL_STORAGE,
                                      parser_->pressure());
  if(reactor_ == NULL) {
    printf("# ERROR: Could not create ConstPressureReactor for files:\n"
           "#            mechanism      file = %s\n"
           "#            thermodynamics file = %s\n"
           "#            log            file = %s\n",
           parser_->mech_file().c_str(),
           parser_->therm_file().c_str(),
           parser_->log_file().c_str());
    exit(-1);
  }
  reactor_->SetReferenceTemperature(parser_->ref_temperature());

  mechanism_ = reactor_->GetMechanism();

  // setup transport interface
  trans_ = transport::InterfaceFactory::CreateMassBased(parser_->transport_model());

  std::vector<std::string> transport_files;
  transport_files.push_back(parser_->trans_file());

  error_code = trans_->Initialize(mechanism_,
                                  transport_files,
                                  parser_->log_file());
  if(error_code != transport::NO_ERROR) {
    printf("# ERROR: Could not Initialize MassTransportInterface for files:\n"
           "#            mechanism      file = %s\n"
           "#            thermodynamics file = %s\n"
           "#            transport      file = %s\n"
           "#            log            file = %s\n"
           "#        Initialize returned error code = %d\n",
           parser_->mech_file().c_str(),
           parser_->therm_file().c_str(),
           parser_->trans_file().c_str(),
           parser_->log_file().c_str(),
           error_code);
    exit(-1);
  }
  // setup logger
  logger_ = new zerork::utilities::Logger(parser_->log_file());
  if(logger_ == NULL) {
    printf("# ERROR: Could not create new logger for file %s\n",
           parser_->log_file().c_str());
    exit(-1);
  }

  SetInlet();
  SetMemory();
  logger_->FFlush();
}

FlameParams::~FlameParams()
{
  if(parser_ != NULL) {
    delete parser_;
  }
  if(reactor_ != NULL) {
    delete reactor_;
  }
  if(trans_ != NULL) {
    delete trans_;
  }
  if(logger_ != NULL) {
    delete logger_;
  }
  if(transport_input_.mass_fraction_ != NULL) {
    delete [] transport_input_.mass_fraction_;
  }
  if(transport_input_.grad_temperature_ != NULL) {
    delete [] transport_input_.grad_temperature_;
  }
  if(transport_input_.grad_pressure_ != NULL) {
    delete [] transport_input_.grad_pressure_;
  }
  if(transport_input_.grad_mass_fraction_ != NULL) {
    delete [] transport_input_.grad_mass_fraction_;
  }
}

// Uses parser_->inlet_fuel_comp(),
//      parser_->inlet_oxidizer_comp(),
//      parser_->inlet_phi() to set inlet_mass_fractions_
void FlameParams::SetInlet()
{
  const int num_species = reactor_->GetNumSpecies();
  std::map<std::string, double>::const_iterator iter;
  std::vector<double> inlet_mole_fractions;
  std::vector<double> fuel_mole_fractions;
  std::vector<double> oxidizer_mole_fractions;
  std::vector<double> exhaust_mole_fractions;
  std::vector<double> exhaust_mass_fractions;
  std::vector<double> molecular_mass;
  double fuel_fraction;
  double phi_term;
  double mole_fraction_sum;

  fuel_mole_fractions.assign(num_species, 0.0);
  oxidizer_mole_fractions.assign(num_species, 0.0);
  inlet_mole_fractions.assign(num_species, 0.0);
  molecular_mass.assign(num_species, 0.0);
  inlet_mass_fractions_.assign(num_species, 0.0);

  exhaust_mole_fractions.assign(num_species, 0.0);
  exhaust_mass_fractions.assign(num_species, 0.0 );

  // fuel compositions
  if(parser_->inlet_fuel_comp().size() == 0) {
    printf("# ERROR: inlet_fuel_comp().size() is zero.\n");
    exit(-1);
  }

  logger_->PrintF("# From inlet_fuel_comp map in %s\n",
                  input_name_.c_str());
  logger_->PrintF("#   Species Name             Mole Fraction\n");

  mole_fraction_sum = 0.0;
  for(iter=parser_->inlet_fuel_comp().begin();
      iter != parser_->inlet_fuel_comp().end();
      ++iter) {

    std::string species_name=iter->first;
    std::string state_name = std::string("MassFraction_")+species_name;
    double mole_fraction = iter->second;
    int species_id = reactor_->GetIdOfState(state_name.c_str());
    logger_->PrintF("%16s  %24.18e\n",
                    species_name.c_str(),
                    mole_fraction);

    if(0 <= species_id  && species_id < num_species) {
      fuel_species_id_.push_back(species_id);
      fuel_mole_fractions[species_id] = mole_fraction;
      mole_fraction_sum += mole_fraction;
    } else {
      printf("# ERROR: did not find species %s in the mechanism\n",
             species_name.c_str());
      exit(-1);
    }
  }
  // re-normalize fuel mole fractions
  NormalizeComposition(num_species,&fuel_mole_fractions[0]);

  // oxidizer compositions
  if(parser_->inlet_oxidizer_comp().size() == 0) {
    printf("# ERROR: inlet_oxidizer_comp().size() is zero.\n");
    exit(-1);
  }

  logger_->PrintF("# From inlet_oxidizer_comp map in %s\n",
                  input_name_.c_str());
  logger_->PrintF("#   Species Name             Mole Fraction\n");

  mole_fraction_sum = 0.0;
  for(iter=parser_->inlet_oxidizer_comp().begin();
      iter != parser_->inlet_oxidizer_comp().end();
      ++iter) {

    std::string species_name=iter->first;
    std::string state_name = std::string("MassFraction_")+species_name;
    double mole_fraction = iter->second;
    int species_id = reactor_->GetIdOfState(state_name.c_str());
    logger_->PrintF("%16s  %24.18e\n",
                    species_name.c_str(),
                    mole_fraction);

    if(0 <= species_id  && species_id <num_species) {
      oxidizer_mole_fractions[species_id] = mole_fraction;
      mole_fraction_sum += mole_fraction;
    } else {
      printf("# ERROR: did not find species %s in the mechanism\n",
             species_name.c_str());
      exit(-1);
    }
  }
  // re-normalize oxidizer mole fractions
  NormalizeComposition(num_species,&oxidizer_mole_fractions[0]);

  // compute the excess atomic oxygen of the fuel and oxidizer stream
  double fuel_atomic_oxygen_sum, oxidizer_atomic_oxygen_sum;
  double excess_atomic_oxygen;
  std::vector<int> num_hydrogen, num_carbon, num_oxygen;
  num_hydrogen.assign(num_species, 0);
  num_carbon.assign(num_species, 0);
  num_oxygen.assign(num_species, 0);

  reactor_->GetSpeciesHydrogenCount(&num_hydrogen[0]);
  reactor_->GetSpeciesCarbonCount(&num_carbon[0]);
  reactor_->GetSpeciesOxygenCount(&num_oxygen[0]);

  fuel_atomic_oxygen_sum = oxidizer_atomic_oxygen_sum = 0.0;
  for(int j=0; j<num_species; ++j) {

    excess_atomic_oxygen =
      -2.0*num_carbon[j] - 0.5*num_hydrogen[j] + 1.0*num_oxygen[j];

    fuel_atomic_oxygen_sum     +=
      fuel_mole_fractions[j]*excess_atomic_oxygen;

    oxidizer_atomic_oxygen_sum +=
      oxidizer_mole_fractions[j]*excess_atomic_oxygen;
  }

  logger_->PrintF(
    "# Moles of excess atomic oxygen/mole of fuel    :  %24.18e\n",
    fuel_atomic_oxygen_sum);
  logger_->PrintF(
    "# Moles of excess atomic oxygen/mole of oxidizer:  %24.18e\n",
    oxidizer_atomic_oxygen_sum);
  logger_->PrintF(
    "# Inlet equivalence ratio                       : %24.18e\n",
    parser_->inlet_phi());

  if(fabs(fuel_atomic_oxygen_sum-oxidizer_atomic_oxygen_sum) < 1.0e-300) {

    printf("# ERROR: can not balance fuel and oxidizer compositions\n"
           "#        because the excess atomic oxygen is the same\n");
    exit(-1);
  }

  phi_term =
    -oxidizer_atomic_oxygen_sum/fuel_atomic_oxygen_sum*parser_->inlet_phi();

  fuel_fraction = phi_term/(1.0+phi_term);
  if(0.0 > fuel_fraction || 1.0 < fuel_fraction) {

    printf("# ERROR: can not balance fuel and oxidizer compositions\n"
           "#        because the fuel fraction is out-of-bounds:\n"
           "#            fuel     excess atomic oxygen = %24.18e\n"
           "#            oxidizer excess atomic oxygen = %24.18e\n"
           "#            fuel fraction = %24.18e\n",
           fuel_atomic_oxygen_sum,
           oxidizer_atomic_oxygen_sum,
           fuel_fraction);
    exit(-1);
  }

  reactor_->GetSpeciesMolecularWeight(&molecular_mass[0]);

  inlet_molecular_mass_ = 0.0;
  for(int j=0; j<num_species; ++j) {

    inlet_mole_fractions[j] =
      fuel_fraction*fuel_mole_fractions[j] +
      (1.0-fuel_fraction)*oxidizer_mole_fractions[j];

    inlet_molecular_mass_ += molecular_mass[j]*inlet_mole_fractions[j];
  }
  // renormalize inlet mole fractions
  NormalizeComposition(num_species,&inlet_mole_fractions[0]);

  for(int j=0; j<num_species; ++j) {
    inlet_mass_fractions_[j] =
      inlet_mole_fractions[j]*molecular_mass[j]/inlet_molecular_mass_;
  }
  // renormalize inlet mass fractions
  NormalizeComposition(num_species,&inlet_mass_fractions_[0]);

  // Add EGR
  double egr = parser_->egr();
  if(egr > 0.0) {
    mechanism_->getMolarIdealExhaust(&inlet_mole_fractions[0],&exhaust_mole_fractions[0]);
    mechanism_->getYfromX(&exhaust_mole_fractions[0],&exhaust_mass_fractions[0]);
    for(int j=0; j<num_species; j++) {
      inlet_mass_fractions_[j] = (1.0-egr)*inlet_mass_fractions_[j]+egr*exhaust_mass_fractions[j];
    }
    mechanism_->getXfromY(&inlet_mass_fractions_[0],&inlet_mole_fractions[0]);

    // Recompute molecular mass
    inlet_molecular_mass_ = 0.0;
    for(int j=0; j<num_species; ++j) {
      inlet_molecular_mass_ += molecular_mass[j]*inlet_mole_fractions[j];
    }
  }

  //if there's inlet_full_comp use that instea
  // fuel compositions
  if(parser_->inlet_full_comp().size() != 0) {
    printf("# Using inlet_full_comp.\n");

    logger_->PrintF("# From inlet_full_comp map in %s\n",
                    input_name_.c_str());
    logger_->PrintF("#   Species Name             Mole Fraction\n");

    mole_fraction_sum = 0.0;
    for(iter=parser_->inlet_full_comp().begin();
        iter != parser_->inlet_full_comp().end();
        ++iter) {

      std::string species_name=iter->first;
      std::string state_name = std::string("MassFraction_")+species_name;
      double mole_fraction = iter->second;
      int species_id = reactor_->GetIdOfState(state_name.c_str());
      logger_->PrintF("%16s  %24.18e\n",
                      species_name.c_str(),
                      mole_fraction);

      if(0 <= species_id  && species_id < num_species) {
        full_species_id_.push_back(species_id);
        inlet_mole_fractions[species_id] = mole_fraction;
        mole_fraction_sum += mole_fraction;
      } else {
        printf("# ERROR: did not find species %s in the mechanism\n",
               species_name.c_str());
        exit(-1);
      }
    }

    // renormalize inlet mole fractions
    NormalizeComposition(num_species,&inlet_mole_fractions[0]);

    // Compute inlet mass fractions
    for(int j=0; j<num_species; ++j) {
      inlet_mass_fractions_[j] =
        inlet_mole_fractions[j]*molecular_mass[j]/inlet_molecular_mass_;
    }
    // renormalize inlet mass fractions
    NormalizeComposition(num_species,&inlet_mass_fractions_[0]);
  }


  logger_->PrintF("# Inlet mixture molecular mass [kg/kmol] = %24.18e\n",
                  inlet_molecular_mass_);
  logger_->PrintF("#   Species Name             Mole Fraction             Mass Fraction\n");

  for(int j=0; j<num_species; ++j) {

    if(inlet_mass_fractions_[j] > 0.0) {
      std::string state_name = reactor_->GetNameOfStateId(j);
      std::string species_name;
      // skip the prefix "MassFraction_" in the state name
      species_name.assign(state_name,13,std::string::npos);
      logger_->PrintF("%16s  %24.18e  %24.18e\n",
                      species_name.c_str(),
                      inlet_mole_fractions[j],
                      inlet_mass_fractions_[j]);
    }
  }
  inlet_temperature_ = parser_->inlet_temperature()/
    parser_->ref_temperature();
  ref_temperature_ =  parser_->ref_temperature();
  pressure_ = parser_->pressure();
  inlet_relative_volume_ =
    reactor_->GetGasConstant()*parser_->inlet_temperature()/
    (parser_->pressure()*inlet_molecular_mass_);
} // void FlameParams::SetInlet()

// requires SetGrid to be set first
void FlameParams::SetMemory()
{
  const int num_species  = reactor_->GetNumSpecies();
  const int num_points = 1; //use a single point to evaluate Lewis numbers

  // Setup the MassTransportFactory data structure
  transport_input_.num_dimensions_        = 1;
  transport_input_.ld_grad_temperature_   = 1;
  transport_input_.ld_grad_pressure_      = 1;
  transport_input_.ld_grad_mass_fraction_ = num_species;
  transport_input_.mass_fraction_         = new double[num_species];
  transport_input_.grad_temperature_      = new double[1]; // num_dimensions
  transport_input_.grad_pressure_         = new double[1]; // num_dimensions
  transport_input_.grad_mass_fraction_    = new double[num_species];

  // Apply constant pressure approximation
  transport_input_.pressure_         = parser_->pressure();
  transport_input_.grad_pressure_[0] = 0.0;

  // create and set the inverse molecular mass array
  inv_molecular_mass_.assign(num_species, 0.0);
  reactor_->GetSpeciesMolecularWeight(&inv_molecular_mass_[0]);
  for(int j=0; j<num_species; ++j) {
    inv_molecular_mass_[j] = 1.0/inv_molecular_mass_[j];
  }

  // create the workspace for the species specific heats
  species_specific_heats_.assign(num_species*(num_points), 0.0);

  // create the workspace for the species mass fluxes and Lewis numbers
  species_mass_flux_.assign(num_species*(num_points+1), 0.0); //larger size for derivatives
  species_lewis_numbers_.assign((num_points+1)*num_species, 0.0);

  // create the workspace for the flux interface conductivities
  thermal_conductivity_.assign(num_points+1, 0.0);//larger size for derivatives

  // create the workspace for the mixture specific heat at each grid point
  mixture_specific_heat_.assign(num_points, 0.0);

  // Are we computing thermodynamic equilibrium
  use_equilibrium_ = parser_->use_equilibrium();

}

static double NormalizeComposition(const size_t num_elements,
                                   double composition[])
{
  double sum = 0.0;
  for(size_t j=0; j<num_elements; ++j) {
    sum += composition[j];
  }
  if(fabs(sum) > 1.0e-300) {
    sum = 1.0/sum;
    for(size_t j=0; j<num_elements; ++j) {
      composition[j] *= sum;
    }
    sum = 1.0/sum;
  }
  return sum;
}
