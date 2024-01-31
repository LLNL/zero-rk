#include <stdlib.h>
#include <stdio.h>
#ifdef _WIN32
#define _USE_MATH_DEFINES //for M_PI
#endif
#include <cmath>

#include <fstream>
#include <utilities/string_utilities.h>
#include <utilities/math_utilities.h>
#include <utilities/file_utilities.h>

#include "flame_params.h"


static double NormalizeComposition(const size_t num_elements,
                                   double composition[]);

FlameParams::FlameParams(const std::string &input_name)
{
  npes_ = 1;
  my_pe_ = 0;
#ifdef ZERORK_MPI
  comm_ = MPI_COMM_WORLD;
  MPI_Comm_size(comm_, &npes_);
  MPI_Comm_rank(comm_, &my_pe_);
#endif

  int error_code;

  parser_  = NULL;
  reactor_ = NULL;
  trans_   = NULL;
  logger_  = NULL;
  sparse_matrix_.clear();
  valid_jacobian_structure_ = true;

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

  // Set inlet values
  SetInlet();
  // Set initial composition
  SetInitialComposition();
  // Set grid
  SetGrid();
  // Set wall temperature
  SetWallProperties();
  // Set arrays
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
  for(size_t j=0; j<sparse_matrix_.size(); ++j) {
    if(sparse_matrix_[j] != NULL) {
      delete sparse_matrix_[j];
    }
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
  std::vector<double> exhaust_mole_fractions;
  std::vector<double> exhaust_mass_fractions;
  std::vector<double> fuel_mole_fractions;
  std::vector<double> oxidizer_mole_fractions;
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
  exhaust_mass_fractions.assign(num_species, 0.0);

  // Get fuel composition
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

  // Get oxidizer composition
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
      oxidizer_species_id_.push_back(species_id);
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

  // Use equivalence ratio to compute mole fractions
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

  // Compute inlet mass fractions
  for(int j=0; j<num_species; ++j) {
    inlet_mass_fractions_[j] =
      inlet_mole_fractions[j]*molecular_mass[j]/inlet_molecular_mass_;
  }
  // renormalize inlet mass fractions
  NormalizeComposition(num_species,&inlet_mass_fractions_[0]);

  // Add EGR to inlet composition
  double egr = parser_->egr();
  if(egr > 0.0) {
    if(parser_->egr_comp().size() == 0) {
      //use ideal EGR formulation if egr_comp is not provided
      mechanism_->getMolarIdealExhaust(&inlet_mole_fractions[0],&exhaust_mole_fractions[0]);
    } else { //
      // Get EGR composition from egr_comp field if provided
      for(int j=0; j<num_species; ++j){exhaust_mole_fractions[j] = 0.0;}

      mole_fraction_sum = 0.0;
      for(iter=parser_->egr_comp().begin();
          iter != parser_->egr_comp().end();
          ++iter) {

        std::string species_name=iter->first;
        std::string state_name = std::string("MassFraction_")+species_name;
        double mole_fraction = iter->second;
        int species_id = reactor_->GetIdOfState(state_name.c_str());

        if(0 <= species_id  && species_id < num_species) {
          egr_species_id_.push_back(species_id);
          exhaust_mole_fractions[species_id] = mole_fraction;
          mole_fraction_sum += mole_fraction;
        } else {
          printf("# ERROR: did not find species %s in the mechanism\n",
                 species_name.c_str());
          exit(-1);
        }
      }

      // renormalize inlet mole fractions
      NormalizeComposition(num_species,&exhaust_mole_fractions[0]);
    }
    mechanism_->getYfromX(&exhaust_mole_fractions[0],&exhaust_mass_fractions[0]);

    for(int j=0; j<num_species; j++) {
      inlet_mass_fractions_[j] = (1.0-egr)*inlet_mass_fractions_[j]+egr*exhaust_mass_fractions[j];
    }
    mechanism_->getXfromY(&inlet_mass_fractions_[0],&inlet_mole_fractions[0]);

    // Recompute
    inlet_molecular_mass_ = 0.0;
    for(int j=0; j<num_species; ++j) {
      inlet_molecular_mass_ += molecular_mass[j]*inlet_mole_fractions[j];
    }
  }

  // If inlet_full_comp is speficied in input file, use that instead
  // fuel compositions
  if(parser_->inlet_full_comp().size() != 0) {
    printf("# Using inlet_full_comp.\n");

    for(int j=0; j<num_species; ++j) {
      inlet_mole_fractions[j] = 0.0;
    }

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

    // Recompute inlet molecular mass
    inlet_molecular_mass_ = 0.0;
    for(int j=0; j<num_species; ++j) {
      inlet_molecular_mass_ += molecular_mass[j]*inlet_mole_fractions[j];
    }

    // Compute inlet mass fractions
    for(int j=0; j<num_species; ++j) {
      inlet_mass_fractions_[j] =
        inlet_mole_fractions[j]*molecular_mass[j]/inlet_molecular_mass_;
    }
    // renormalize inlet mass fractions
    NormalizeComposition(num_species,&inlet_mass_fractions_[0]);

  }//if inlet_full_comp

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

  inlet_relative_volume_ =
    reactor_->GetGasConstant()*parser_->inlet_temperature()/
    (parser_->pressure()*inlet_molecular_mass_);

} // void FlameParams::SetInlet()

// Set initial composition vector if initial_comp is defined in input file
void FlameParams::SetInitialComposition()
{
  const int num_species = reactor_->GetNumSpecies();
  std::map<std::string, double>::const_iterator iter;
  std::vector<double> initial_mole_fractions;
  std::vector<double> molecular_mass;
  double mass_fraction_sum = 0.0;
  double mixture_molecular_mass = 0.0;

  initial_mass_fractions_.assign(num_species, 0.0);
  initial_mole_fractions.assign(num_species, 0.0);
  molecular_mass.assign(num_species, 0.0);

  reactor_->GetSpeciesMolecularWeight(&molecular_mass[0]);

  // initial compositions
  if(parser_->initial_comp().size() == 0) {
    printf("# ERROR: initial_comp().size() is zero.\n");
    exit(-1);
  }

  logger_->PrintF("# From initial_comp map in %s\n",
                  input_name_.c_str());
  logger_->PrintF("#   Species Name             Mole Fraction\n");

  mixture_molecular_mass = 0.0;
  for(iter=parser_->initial_comp().begin();
      iter != parser_->initial_comp().end();
      ++iter) {

    std::string species_name=iter->first;
    std::string state_name = std::string("MassFraction_")+species_name;
    double mole_fraction = iter->second;
    int species_id = reactor_->GetIdOfState(state_name.c_str());

    if(0 <= species_id  && species_id <num_species) {
      initial_mole_fractions[species_id] = mole_fraction;
      // compute the mixture molecular mass (to within a scaling factor
      // equal to the mole fraction sum)
      mixture_molecular_mass += mole_fraction*molecular_mass[species_id];

      logger_->PrintF("%16s  %24.18e\n",
                      species_name.c_str(),
                      mole_fraction);
    } else {
      printf("# ERROR: In FlameParams::SetInitialComposition,\n"
             "#        did not find species %s given by initial_comp\n"
             "#        in the mechanism\n",
             species_name.c_str());
      exit(-1);
    }
  }

  mass_fraction_sum = 0.0;
  for(int j=0; j<num_species; ++j) {
    initial_mass_fractions_[j] = initial_mole_fractions[j]*
      (molecular_mass[j]/mixture_molecular_mass);
    mass_fraction_sum += initial_mass_fractions_[j];
  }

  // normalize the initial mass fractions
  if(mass_fraction_sum > 1.0e-300) {
    for(int j=0; j<num_species; ++j) {
      initial_mass_fractions_[j] /= mass_fraction_sum;
    }
  }
}

// Set grid
void FlameParams::SetGrid()
{
  length_ = parser_->length();
  int num_points = parser_->num_points();

  if(parser_->grid_file() == std::string(zerork::utilities::null_filename)) {
    // if no grid provided -> uniform grid from input file parameters
    if(parser_->num_points() < 2) {
      printf("# ERROR: number of grid points must be two or greater than two\n"
	     "#            num_points = %d\n", parser_->num_points());
      exit(-1);
    }
    if(parser_->length() <= 0.0) {
      printf("# ERROR: flame tube length must be positive\n"
	     "#            length = %.18g\n",
	     parser_->length());
      exit(-1);
    }

    const double delta_z = length_/(double)(num_points+1-1);
    num_points_ = num_points;
    z_.assign(num_points, 0.0);
    for(int j=0; j<num_points; ++j) {
      z_[j] = delta_z + (double)j*delta_z;
    }

  } else {
    // Otherwise, read grid from file
    // Right now all processors are doing the same read here
    std::ifstream grid_file;
    std::string line;
    std::string delimiters = std::string(zerork::utilities::WHITESPACE) + ",";
    std::vector<std::string> fields;

    grid_file.open(parser_->grid_file().c_str());

    if(grid_file) {
      while(zerork::utilities::GetAnyLine(grid_file,&line)) {
	zerork::utilities::SplitStringToVector(line,
                                                delimiters,
                                                &fields);
        if(fields.size() == 1) {
          if(zerork::utilities::StringIsDouble(fields[0])) {
            z_.push_back(atof(fields[0].c_str()));
          }
        }
      }

    } else {
      printf("# ERROR: could not open grid file = %s\n",
             parser_->grid_file().c_str());
      exit(-1);
    }//if file

    if(z_.size() == 0) {
      printf("# ERROR: no position data found in\n"
             "#        grid file = %s\n",
             parser_->grid_file().c_str());
      exit(-1);
    } // if z.size

    length_ = z_[z_.size()-1];
    //Remove first point
    z_.erase(z_.begin());
  }

  num_points = z_.size();
  num_points_ = num_points;
  num_local_points_ = num_points/npes_;

  if(num_points % npes_ != 0 ) {
    printf("Number of grid points not divisible by number of processors \n");
#ifdef ZERORK_MPI
    MPI_Finalize();
#endif
    exit(-1);
  }

  // Compute midpoints (zm) and spacings (dz, dzm)
  zm_.assign(num_points, 0.0);
  dz_.assign(num_points, 1.0);
  dzm_.assign(num_points, 1.0);

  int nover = 2;
  if(num_local_points_ < nover ) {
    printf("Need at least two grid points per processor for second order discretization \n");
#ifdef ZERORK_MPI
    MPI_Finalize();
#endif
    exit(-1);
  }

#ifdef ZERORK_MPI
  MPI_Status status;
#endif
  dz_local_.assign( num_local_points_+(2*nover), 0.0);
  dzm_local_.assign( num_local_points_+(2*nover), 0.0);
  inv_dz_local_.assign( num_local_points_+(2*nover), 0.0);
  inv_dzm_local_.assign( num_local_points_+(2*nover), 0.0);

  // Compute midpoints and grid spacings
  for(int j=1; j<num_points; ++j) {
    dz_[j] = z_[j]-z_[j-1];
    zm_[j] = 0.5*(z_[j]+z_[j-1]);
  }
  dz_[0] = dz_[1]; //hack
  zm_[0] = z_[0]-0.5*dz_[0];
  for(int j=0; j<num_points-1; ++j) {
    dzm_[j] = zm_[j+1]-zm_[j];
  }
  dzm_[num_points-1] = dz_[num_points-1];

  for (int j=0; j<num_local_points_; ++j) {
    int jglobal = j + my_pe_*num_local_points_;
    dz_local_[nover+j] = dz_[jglobal];
    dzm_local_[nover+j] = dzm_[jglobal];
    inv_dz_local_[nover+j] = 1.0/dz_[jglobal];
    inv_dzm_local_[nover+j] = 1.0/dzm_[jglobal];
  }

#ifdef ZERORK_MPI
  // Send left and right
  if (my_pe_ !=0) {
    MPI_Send(&dz_local_[nover], nover, MPI_DOUBLE, my_pe_-1, 0, comm_);
    MPI_Send(&dzm_local_[nover], nover, MPI_DOUBLE, my_pe_-1, 0, comm_);
    MPI_Send(&inv_dz_local_[nover], nover, MPI_DOUBLE, my_pe_-1, 0, comm_);
    MPI_Send(&inv_dzm_local_[nover], nover, MPI_DOUBLE, my_pe_-1, 0, comm_);
  }
  if (my_pe_ != npes_-1) {
    MPI_Send(&dz_local_[num_local_points_], nover, MPI_DOUBLE, my_pe_+1, 0, comm_);
    MPI_Send(&dzm_local_[num_local_points_], nover, MPI_DOUBLE, my_pe_+1, 0, comm_);
    MPI_Send(&inv_dz_local_[num_local_points_], nover, MPI_DOUBLE, my_pe_+1, 0, comm_);
    MPI_Send(&inv_dzm_local_[num_local_points_], nover, MPI_DOUBLE, my_pe_+1, 0, comm_);
  }

  if (my_pe_ !=0) {
    MPI_Recv(&dz_local_[0], nover, MPI_DOUBLE, my_pe_-1, 0, comm_, &status);
    MPI_Recv(&dzm_local_[0], nover, MPI_DOUBLE, my_pe_-1, 0, comm_, &status);
    MPI_Recv(&inv_dz_local_[0], nover, MPI_DOUBLE, my_pe_-1, 0, comm_, &status);
    MPI_Recv(&inv_dzm_local_[0], nover, MPI_DOUBLE, my_pe_-1, 0, comm_, &status);
  }
  if (my_pe_ != npes_-1) {
    MPI_Recv(&dz_local_[num_local_points_+nover], nover, MPI_DOUBLE, my_pe_+1, 0, comm_, &status);
    MPI_Recv(&dzm_local_[num_local_points_+nover], nover, MPI_DOUBLE, my_pe_+1, 0, comm_, &status);
    MPI_Recv(&inv_dz_local_[num_local_points_+nover], nover, MPI_DOUBLE, my_pe_+1, 0, comm_, &status);
    MPI_Recv(&inv_dzm_local_[num_local_points_+nover], nover, MPI_DOUBLE, my_pe_+1, 0, comm_, &status);
  }
#endif

  if(my_pe_ == 0) {
    for(int j=0; j<nover; ++j) {
      dz_local_[j] = dz_[0];
      dzm_local_[j] = dzm_[0];
      inv_dz_local_[j] = 1.0/dz_[0];
      inv_dzm_local_[j] = 1.0/dzm_[0];
    }
  }

  if(my_pe_ == npes_-1) {
    for(int j=num_local_points_+nover; j<num_local_points_+2*nover; ++j) {
      dz_local_[j] = dz_[num_points-1];
      dzm_local_[j] = dzm_[num_points-1];
      inv_dz_local_[j] = 1.0/dz_[num_points-1];
      inv_dzm_local_[j] = 1.0/dzm_[num_points-1];
    }
  }

  // set mass flux
  inlet_relative_volume_ =
    reactor_->GetGasConstant()*parser_->inlet_temperature()/
    (parser_->pressure()*inlet_molecular_mass_);

  mass_flux_.assign(num_local_points_,0.0);
  // for flame in tube simulations
  if(parser_->simulation_type() == 0) {
    mass_flux_inlet_ = 4.0*parser_->mass_flow()/
      (M_PI*parser_->diameter()*parser_->diameter());
  } else { // for flame speed calculations
    mass_flux_inlet_ = parser_->inlet_velocity()/
      (reactor_->GetGasConstant()*parser_->inlet_temperature()/
       (parser_->pressure()*inlet_molecular_mass_));
  }

  // Initialize mass flux to inlet value
  for (int j=0; j<num_local_points_; ++j) {
    mass_flux_[j] = mass_flux_inlet_;
  }

  diameter_ = parser_->diameter();
}

// Set wall properties for flame in tube simulations
void FlameParams::SetWallProperties()
{
  const int num_points = z_.size();
  std::vector<double> file_wall_position, file_wall_temperature;

  wall_temperature_.clear();
  wall_temperature_.assign(num_points,
                           parser_->inlet_temperature()/
                           parser_->ref_temperature());

  if(parser_->wall_temperature_file() == std::string(zerork::utilities::null_filename)) {
    // if no wall temperature file is specified, then wall temperature
    // remains equal to the inlet_temperature and the nusselt number is
    // set to zero for an adiabatic calculation
    nusselt_ = 0.0;
  } else {
    nusselt_ = parser_->nusselt();

    std::ifstream wall_temperature_file;
    std::string line;
    std::string delimiters = std::string(zerork::utilities::WHITESPACE) + ",";
    std::vector<std::string> fields;
    wall_temperature_file.open(parser_->wall_temperature_file().c_str());

    if(wall_temperature_file) {
      while(zerork::utilities::GetAnyLine(wall_temperature_file,&line)) {
        // remove everything before the comment character
	zerork::utilities::SplitBeforeDelimiter(line,
                                        "!",
                                        true,
                                        &line);
	zerork::utilities::SplitStringToVector(line,
                                       delimiters,
                                       &fields);
        // only record valid data lines
        if(fields.size() == 2) {
          if(zerork::utilities::StringIsDouble(fields[0]) &&
             zerork::utilities::StringIsDouble(fields[1])) {

            file_wall_position.push_back(atof(fields[0].c_str()));
            file_wall_temperature.push_back(atof(fields[1].c_str()));
          }
        }
      } // while(zerork::utilities::GetAnyLine(wall_temperature_file,&line))

    } else {
      printf("# ERROR: could not open wall_temperature_file = %s\n",
             parser_->wall_temperature_file().c_str());
      exit(-1);
    }

    if(file_wall_position.size() == 0) {
      printf("# ERROR: no wall position and temperature data found in\n"
             "#        wall_temperature_file = %s\n",
             parser_->wall_temperature_file().c_str());
      exit(-1);
    }

    zerork::utilities::SortVectors(&file_wall_position, &file_wall_temperature);
    logger_->PrintF("# Wall temperature profile from file = %s\n",
                    parser_->wall_temperature_file().c_str());
    logger_->PrintF("# Number of file points = %d\n",
                    (int)file_wall_position.size());
    logger_->PrintF("# Position [m]              Wall Temperature [K]\n");
    for(size_t j=0; j<file_wall_position.size(); ++j) {
      logger_->PrintF("%24.16e  %24.16e\n",
                      file_wall_position[j],
                      file_wall_temperature[j]);
    }


    logger_->PrintF("# Wall temperature profile on the computational grid\n");
    logger_->PrintF("# Number of grid points = %d\n",(int)z_.size());
    logger_->PrintF("# Position [m]              Wall Temperature [K]\n");
    zerork::utilities::InterpolationTable table(file_wall_position.size(),
                                        &file_wall_position[0],
                                        &file_wall_temperature[0],
                                        zerork::utilities::LINEAR,
                                        false); // use_interpolation
    const int num_points = z_.size();
    for(int j=0; j<num_points; ++j) {
      wall_temperature_[j] = table.Interpolate(z_[j]);
      logger_->PrintF("%24.16e  %24.16e\n",
                      z_[j],
                      wall_temperature_[j]);
      wall_temperature_[j] /= parser_->ref_temperature();
    }

  } // if(parser_->wall_temperature_file() == std::string(zerork::utilities::null_filename)) else

}

// requires SetGrid to be set first
void FlameParams::SetMemory()
{
  const int num_species  = reactor_->GetNumSpecies();
  const int num_states   = reactor_->GetNumStates();
  const int num_nonzeros = reactor_->GetJacobianSize();
  const int num_reactors = num_local_points_;
  const int num_local_points = num_local_points_;

  std::vector<int> column_id;
  int last_dense_id, dense_id;

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

  // create extended (with ghost cells) state variable array
  // for parallel simulations. nover=2 for second order operators
  int nover = 2;
  y_ext_.assign( (num_local_points_+(2*nover))*num_states, 0.0);
  mass_flux_ext_.assign( num_local_points_+(2*nover), 0.0);

  // create the workspace for the species specific heats
  species_specific_heats_.assign(num_species*(num_local_points+1), 0.0);

  // create the workspace for the species mass fluxes and Lewis numbers
  species_mass_flux_.assign(num_species*(num_local_points+1), 0.0); //larger size for derivatives
  species_lewis_numbers_.assign((num_local_points+1)*num_species, 0.0);

  // create the workspace for the flux interface conductivities
  thermal_conductivity_.assign(num_local_points+1, 0.0);//larger size for derivatives

  // create the workspace for the mixture specific heat at each grid point
  mixture_specific_heat_.assign(num_local_points+1, 0.0);

  // create the workspace for the mixture specific heat at each interface
  mixture_specific_heat_mid_.assign(num_local_points+1, 0.0);

  // Get convective scheme type
  convective_scheme_type_ = parser_->convective_scheme_type();

  // Get simulation type (0 for tube, 1 for flame speed)
  simulation_type_ = parser_->simulation_type();

  // create the Jacobian data structures needed for the banded and SPGMR
  // integrators
  store_jacobian_  = parser_->store_jacobian();

  // Reaction rate limiter
  step_limiter_.assign( reactor_->GetNumSteps(), parser_->step_limiter() );

  // For tracking/output
  max_velocity_ = 0.0;
  min_velocity_ = 0.0;
  max_thermal_diffusivity_ = 0.0;

  // Setup sparse Jacobian arrays
  row_id_.assign(num_nonzeros, 0);
  column_id.assign(num_nonzeros, 0);
  column_sum_.assign(num_states+1,0);

  reactor_jacobian_.assign(num_nonzeros, 0.0);

  // The Jacobian pattern is assumed to be in compressed column storage
  reactor_->GetJacobianPattern(&row_id_[0],
                               &column_id[0]);

  last_dense_id = -1;
  for(int j=0; j<num_nonzeros; ++j) {

    dense_id = num_states*column_id[j] + row_id_[j];

    // check that the column index and row index are increasing consistent
    // with the compressed column storage pattern
    if(dense_id > last_dense_id) {

      last_dense_id = dense_id;
      ++column_sum_[column_id[j]+1];

    } else {

      logger_->PrintF(
        "# ERROR: reactor Jacobian element %d is not in compressed column"
        "#        storage order:\n"
        "#              J element %d (row, col) = (%d, %d)\n"
        "#              J element %d (row, col) = (%d, %d)\n"
        "#        Jacobian structure invalid.\n",
        j,
        j-1,
        row_id_[j-1],
        column_id[j-1],
        j,
        row_id_[j],
        column_id[j],
        j);

      valid_jacobian_structure_ = false;
      exit(-1); // TODO: add recoverable failure

    }  // if(dense_id > last_dense_id) else

  } // for(int j=0; j<num_nonzeros; ++j)

    // At this point the column sum contains just the count per column, need
    // to loop through to change the count to a sum
  for(int j=0; j<num_states; ++j) {
    column_sum_[j+1] += column_sum_[j];
  }

  diagonal_id_.assign(num_states, -1); // use negative one as flag to
                                       // check if there are any missing
                                       // diagonal terms that would prevent
                                       // its use for the sparse factorization
  for(int j=0; j<num_nonzeros; ++j) {
    if(row_id_[j] == column_id[j]) {
      diagonal_id_[row_id_[j]] = j;
    }
  }
  for(int j=0; j<num_states; ++j) {
    if(diagonal_id_[j] == -1) {
      logger_->PrintF(
        "# ERROR: reactor Jacobian does not have diagonal in row %d\n",
        j);

      valid_jacobian_structure_ = false;
      exit(-1); // TODO: add recoverable failure
    }
  }

  // Create the vector of the sparse matrix elements
  sparse_matrix_.assign(num_reactors, NULL);
  for(int j=0; j<num_reactors; ++j) {
    sparse_matrix_[j] = new SparseMatrix(num_states, num_nonzeros);
    if(sparse_matrix_[j] == NULL) {
      logger_->PrintF(
        "# ERROR: failed to allocate new SparseMatrix object\n"
        "#        for grid point %d (z = %.18g [m])\n",
          j, z_[j]);
      exit(-1); // TODO: add recoverable failure
    }
  }

  if(store_jacobian_) {
    // TODO: add allocation check in case we run out of memory for large
    //       mechanisms and grids
    saved_jacobian_.assign(num_nonzeros*num_reactors, 0.0);
  }

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
