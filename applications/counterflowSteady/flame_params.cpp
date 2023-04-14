#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
#include <utilities/string_utilities.h>
#include <utilities/math_utilities.h>
#include <utilities/file_utilities.h>

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
  sparse_matrix_ = NULL;
  sparse_matrix_dist_ = NULL;
  sparse_matrix_chem_.clear();
  valid_jacobian_structure_ = true;
  valid_jacobian_structure_ = true;

  if(!zerork::utilities::FileIsReadable(input_name)) {
    printf("# ERROR: Input file %s is not readable\n",input_name.c_str());
    exit(-1);
  }

  // setup parser
  input_name_ = input_name;
  parser_ = new SteadyFlameIFP(input_name);
  if(parser_ == NULL) {
    printf("# ERROR: Parser for file %s is not created\n",input_name.c_str());
    exit(-1);
  }

  comm_ = MPI_COMM_WORLD;
  inter_comm_ = MPI_COMM_WORLD;
  MPI_Comm_size(comm_, &npes_);
  MPI_Comm_rank(comm_, &my_pe_);
  comm_rank_ = 0;
  num_comms_ = 1;
  if(parser_->sensitivity_analysis() && 
    parser_->sensitivity_processors_per_solution() > 0) {

    if(npes_ % parser_->sensitivity_processors_per_solution() !=0 ) {
      printf("Warning: Total number of processors not evenly divisible by sensitivity_processors_per_solution\n");
      printf("         Continuing without parallelizing sensitivity analysis\n");
    } else if(npes_ / parser_->sensitivity_processors_per_solution() <= 1) {
      printf("Warning: Total number of processors less than or equal to sensitivity_processors_per_solution\n");
      printf("         Continuing without parallelizing sensitivity analysis\n");
    } else {
      num_comms_ = npes_ / parser_->sensitivity_processors_per_solution();
      comm_rank_ = my_pe_ / parser_->sensitivity_processors_per_solution();
      MPI_Comm_split(MPI_COMM_WORLD, comm_rank_, my_pe_, &comm_);
      MPI_Comm_rank(comm_, &my_pe_);
      MPI_Comm_size(comm_, &npes_);
      int inter_comm_rank = 0;
      if(my_pe_ != 0) inter_comm_rank = 1;
      MPI_Comm_split(MPI_COMM_WORLD, inter_comm_rank, comm_rank_, &inter_comm_);
    }
  }

  // setup constant pressure reactor
  reactor_ = new CounterflowReactor(parser_->mech_file().c_str(),
                                    parser_->therm_file().c_str(),
                                    parser_->log_file().c_str(),
                                    COMPRESSED_COL_STORAGE,
                                    parser_->pressure(),
                                    parser_->finite_separation());
  if(reactor_ == NULL) {
    printf("# ERROR: Could not create CounterflowReactor for files:\n"
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
  // Set grid
  SetGrid();
  // Set arrays
  SetMemory();

  logger_->FFlush();

}

FlameParams::~FlameParams() {
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
  if(integrator_type_ == 2) {
    if (superlu_serial_) {
      if (sparse_matrix_ != NULL) {
        delete sparse_matrix_;
      }
    } else {
      if (sparse_matrix_dist_ != NULL) {
        delete sparse_matrix_dist_;
      }
    }
  }
  if(integrator_type_ == 3) {
    for(size_t j=0; j<sparse_matrix_chem_.size(); ++j) {
      if(sparse_matrix_chem_[j] != NULL) {
        delete sparse_matrix_chem_[j];
      }
    }
  }
  if(num_comms_> 1) {
    MPI_Comm_free(&comm_);
    MPI_Comm_free(&inter_comm_);
  }
}

// Uses parser_->inlet_fuel_comp(),
//      parser_->inlet_oxidizer_comp(),
void FlameParams::SetInlet()
{
  const int num_species = reactor_->GetNumSpecies();
  std::map<std::string, double>::const_iterator iter;
  std::vector<double> exhaust_mole_fractions;
  std::vector<double> exhaust_mass_fractions;
  std::vector<double> fuel_mole_fractions;
  std::vector<double> oxidizer_mole_fractions;
  std::vector<double> stoichiometric_mole_fractions;
  std::vector<double> molecular_mass;
  std::vector<double> inlet_mole_fractions;
  double fuel_fraction;
  double phi_term;
  double mole_fraction_sum;

  fuel_mole_fractions.assign(num_species, 0.0);
  oxidizer_mole_fractions.assign(num_species, 0.0);
  stoichiometric_mole_fractions.assign(num_species, 0.0);
  inlet_mole_fractions.assign(num_species, 0.0);
  molecular_mass.assign(num_species, 0.0);

  fuel_mass_fractions_.assign(num_species, 0.0);
  oxidizer_mass_fractions_.assign(num_species, 0.0);
  stoichiometric_mass_fractions_.assign(num_species, 0.0);
  inlet_mass_fractions_.assign(num_species, 0.0);

  exhaust_mole_fractions.assign(num_species, 0.0);
  exhaust_mass_fractions.assign(num_species, 0.0);

  reactor_->GetSpeciesMolecularWeight(&molecular_mass[0]);

  // Get fuel composition
  if(parser_->inlet_fuel_comp().size() == 0) {
    printf("# ERROR: inlet_fuel_comp().size() is zero.\n");
    exit(-1);
  }

  logger_->PrintF("# From inlet_fuel_comp map in %s\n",
                  input_name_.c_str());
  logger_->PrintF("#   Species Name             Mole Fraction\n");

  mole_fraction_sum = 0.0;
  fuel_molecular_mass_ = 0.0;
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
      fuel_molecular_mass_ += mole_fraction*molecular_mass[species_id];
    } else {
      printf("# ERROR: did not find species %s in the mechanism\n",
             species_name.c_str());
      exit(-1);
    }
  }
  // re-normalize fuel mole fractions
  NormalizeComposition(num_species,&fuel_mole_fractions[0]);

  // Compute fuel mass fractions
  for(int j=0; j<num_species; ++j) {
    fuel_mass_fractions_[j] =
      fuel_mole_fractions[j]*molecular_mass[j]/fuel_molecular_mass_;
  }
  // renormalize fuel mass fractions
  NormalizeComposition(num_species,&fuel_mass_fractions_[0]);

  // Get oxidizer composition
  if(parser_->inlet_oxidizer_comp().size() == 0) {
    printf("# ERROR: inlet_oxidizer_comp().size() is zero.\n");
    exit(-1);
  }

  logger_->PrintF("# From inlet_oxidizer_comp map in %s\n",
                  input_name_.c_str());
  logger_->PrintF("#   Species Name             Mole Fraction\n");

  mole_fraction_sum = 0.0;
  oxidizer_molecular_mass_ = 0.0;
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
      oxidizer_molecular_mass_ += mole_fraction*molecular_mass[species_id];
    } else {
      printf("# ERROR: did not find species %s in the mechanism\n",
             species_name.c_str());
      exit(-1);
    }
  }
  // re-normalize oxidizer mole fractions
  NormalizeComposition(num_species,&oxidizer_mole_fractions[0]);

  // Compute oxidizer mass fractions
  for(int j=0; j<num_species; ++j) {
    oxidizer_mass_fractions_[j] =
      oxidizer_mole_fractions[j]*molecular_mass[j]/oxidizer_molecular_mass_;
  }
  // renormalize oxidizer mass fractions
  NormalizeComposition(num_species,&oxidizer_mass_fractions_[0]);


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


  // Add EGR
  double egr = parser_->egr();
  if(egr > 0.0) {
    mechanism_->getMolarIdealExhaust(&inlet_mole_fractions[0],&exhaust_mole_fractions[0]);
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

  // If inlet_full_comp is speficied in input file, use that instead of
  // fuel+oxidizer compositions
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


  // Compute stoichiometric composition
  stoichiometric_molecular_mass_ = 0.0;
  for(int j=0; j<num_species; ++j) {

    stoichiometric_mole_fractions[j] =
      fuel_fraction*fuel_mole_fractions[j] +
      (1.0-fuel_fraction)*oxidizer_mole_fractions[j];

    stoichiometric_molecular_mass_ += molecular_mass[j]*stoichiometric_mole_fractions[j];
  }
  // renormalize stoich mole fractions
  NormalizeComposition(num_species,&stoichiometric_mole_fractions[0]);

  // Compute stoichiometric mass fractions
  stoichiometric_mixture_fraction_ = 0.0;
  for(int j=0; j<num_species; ++j) {
    stoichiometric_mass_fractions_[j] =
      fuel_fraction*fuel_mole_fractions[j]*molecular_mass[j]/stoichiometric_molecular_mass_;
    stoichiometric_mixture_fraction_ += stoichiometric_mass_fractions_[j];
  }
  // Don't renormalize!

  ref_temperature_ =  parser_->ref_temperature();
  ref_momentum_ = parser_->ref_momentum();

  fuel_temperature_ = parser_->fuel_temperature()/ref_temperature_;
  oxidizer_temperature_ = parser_->oxidizer_temperature()/ref_temperature_;

  fuel_relative_volume_ =
    reactor_->GetGasConstant()*parser_->fuel_temperature()/
    (parser_->pressure()*fuel_molecular_mass_);

  oxidizer_relative_volume_ =
    reactor_->GetGasConstant()*parser_->oxidizer_temperature()/
    (parser_->pressure()*oxidizer_molecular_mass_);

  inlet_relative_volume_ =
    reactor_->GetGasConstant()*parser_->fuel_temperature()/ //inlet_temp = fuel_temp
    (parser_->pressure()*inlet_molecular_mass_);

  strain_rate_ = parser_->strain_rate(); //only used with infinite separation

} // void FlameParams::SetInlet()


// Set grid
void FlameParams::SetGrid()
{
  // Initialize with values from input file, overwrite later if necessary
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
      printf("# ERROR: flame domain length must be positive\n"
	     "#            length = %.18g\n",
	     parser_->length());
      exit(-1);
    }

    const double delta_z = length_/(double)(num_points+2-1);
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
    //Remove first and last points (0 and 1)
    z_.erase(z_.begin());
    z_.pop_back();
  }

  num_points = z_.size();
  num_points_ = num_points;
  num_local_points_ = num_points/npes_;

  if(num_points % npes_ != 0 ) {
    printf("Number of grid points not divisible by number of processors \n");
    MPI_Finalize();
    exit(-1);
  }

  // Compute midpoints (zm) and spacings (dz, dzm)
  zm_.assign(num_points, 0.0);
  dz_.assign(num_points, 1.0);
  dzm_.assign(num_points, 1.0);

  int nover = 2;
  nover_ = nover;
  if(num_local_points_ < nover ) {
    printf("Need at least two grid points per processor for second order discretization \n");
    MPI_Finalize();
    exit(-1);
  }

  MPI_Status status;
  dz_local_.assign( num_local_points_+(2*nover), 0.0);
  dzm_local_.assign( num_local_points_+(2*nover), 0.0);
  inv_dz_local_.assign( num_local_points_+(2*nover), 0.0);
  inv_dzm_local_.assign( num_local_points_+(2*nover), 0.0);

  // Compute midpoints and grid spacings
  for(int j=1; j<num_points; ++j) {
    dz_[j] = z_[j]-z_[j-1];
    zm_[j] = 0.5*(z_[j]+z_[j-1]);
  }
  //dz_[0] = dz_[1]; //hack
  dz_[0] = z_[0];
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
      //dz_local_[j] = dz_[num_points-1];
      dz_local_[j] = length_-z_[num_points-1];
      dzm_local_[j] = dzm_[num_points-1];
      inv_dz_local_[j] = 1.0/dz_[num_points-1];
      inv_dzm_local_[j] = 1.0/dzm_[num_points-1];
    }
  }

}

// requires SetGrid to be set first
void FlameParams::SetMemory()
{
  const int num_species  = reactor_->GetNumSpecies();
  const int num_states   = reactor_->GetNumStates();
  const int num_nonzeros_zerod = reactor_->GetJacobianSize();
  const int num_reactors = num_local_points_;
  const int num_points   = z_.size();
  const int num_local_points = num_local_points_;
  int num_local_states = num_local_points*num_states;
  int num_total_states = num_points*num_states;
  num_nonzeros_loc_ = 0;

  std::vector<int> row_id;
  int last_dense_id, dense_id;
  std::vector<int> row_id_zerod, col_id_zerod;

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
  mass_flux_.assign( num_local_points_, 0.0);
  mass_flux_ext_.assign( num_local_points_+(2*nover), 0.0);
  rhs_ext_.assign( (num_local_points_+(2*nover_))*num_states, 0.0);
  rhsConv_.assign(num_local_points_*num_states,0.0);
  rel_vol_.assign(num_local_points_,0.0);
  rel_vol_ext_.assign(num_local_points_+2*nover_,0.0);

  // "old" state vector for pseudo unsteady
  y_old_.assign(num_local_points_*num_states, 0.0);

  // create the workspace for the species specific heats
  species_specific_heats_.assign(num_species*(num_local_points+1), 0.0);

  // create the workspace for the species mass fluxes and Lewis numbers
  species_mass_flux_.assign(num_species*(num_local_points+1), 0.0); //larger size for derivatives
  species_lewis_numbers_.assign((num_local_points+1)*num_species, 0.0);

  // create the workspace for the flux interface conductivities
  thermal_conductivity_.assign(num_local_points+1, 0.0);//larger size for derivatives

  // create the workspace for the mixture viscosity at each interface
  mixture_viscosity_.assign(num_local_points+1, 0.0);//larger size for derivatives

  // create the workspace for the mixture specific heat at each grid point
  mixture_specific_heat_.assign(num_local_points+1, 0.0);

  // create the workspace for the mixture specific heat at each interface
  mixture_specific_heat_mid_.assign(num_local_points+1, 0.0);

  // create the workspace for the mixture molecular mass at each interface
  molecular_mass_mix_mid_.assign(num_local_points+1, 0.0);

  // create the workspace for the fixed temperature profile
  fixed_temperature_.assign(num_local_points, 0.0);

  // Get convective scheme type
  convective_scheme_type_ = parser_->convective_scheme_type();

  // Get number of off-diagonals terms to keep in block jacobian
  num_off_diagonals_ = nover_*num_states;

  // Get simulation type
  simulation_type_ = parser_->simulation_type();

  // Get flame type
  flame_type_ = parser_->flame_type();

  // create the Jacobian data structures needed for the banded and SPGMR
  // integrators
  store_jacobian_  = parser_->store_jacobian();

  // Reaction rate limiter
  step_limiter_.assign( reactor_->GetNumSteps(), parser_->step_limiter() );

  // Pseudo-unsteady solver
  pseudo_unsteady_ = parser_->pseudo_unsteady();


  // For tracking/output
  max_velocity_ = 0.0;
  min_velocity_ = 0.0;
  max_thermal_diffusivity_ = 0.0;
  max_svf_ = 0.0;

  stagnation_plane_ = length_*0.25; //initialize at a quarter?
  jContBC_ = num_points / 4;

  if(parser_->finite_separation()) {
    if(flame_type_ == 0 || flame_type_ == 2) {
      mass_flux_fuel_ = parser_->mass_flux_fuel();
      mass_flux_oxidizer_ = -parser_->mass_flux_oxidizer();
    } else if (flame_type_ == 1) {
      mass_flux_fuel_ = parser_->mass_flux_fuel();
      mass_flux_oxidizer_ = 0.0;
    }
  } else {
    mass_flux_fuel_ = strain_rate_*stagnation_plane_/fuel_relative_volume_;
    mass_flux_oxidizer_ = -strain_rate_*(length_-stagnation_plane_)/oxidizer_relative_volume_;
  }

  // Use SuperLU serial
  superlu_serial_ = parser_->superlu_serial();
  if(npes_ > 1 && superlu_serial_) {
    superlu_serial_ = false;
    if(my_pe_ == 0)
      printf("Can't use serial SuperLU when np > 1. Using distributed SuperLU.\n");
  }

  // create the Jacobian data structures needed for the banded and SPGMR
  // integrators
  integrator_type_ = parser_->integrator_type();
  store_jacobian_  = parser_->store_jacobian();

  // Set Jacobian parameters
  // Case 1 is a block tridiagonal (potentially sparse) matrix solved with SuperLU
  if(integrator_type_ == 2) {
    row_id_zerod.assign(num_nonzeros_zerod, 0);
    col_id_zerod.assign(num_nonzeros_zerod, 0);

    // The Jacobian pattern is assumed to be in compressed column storage
    reactor_->GetJacobianPattern(&row_id_zerod[0],
                                 &col_id_zerod[0]);

    // Set to 1 to force dense blocks
    dense_to_sparse_.assign(num_states*num_states, 1); //, 1); //force dense!!
    dense_to_sparse_offdiag_.assign(num_states*num_states, 1);
    printf("# WARNING: Using dense block tridiagonal Jacobian\n");
    // The following is to use sparse blocks, only matters if
    // dense_to_sparse_/offdiag_ is initialized at 0 above

    // Chemical jacobian pattern -- only for diagonal block
    for (int j=0; j<num_nonzeros_zerod; j++) {
      dense_id = num_states*col_id_zerod[j] + row_id_zerod[j];
      dense_to_sparse_[dense_id] = 1;
      //dense_to_sparse_offdiag_[dense_id] = 1; //for off-diagonal too?
    }

    int num_fuel_species = fuel_species_id_.size();
    int num_oxidizer_species = oxidizer_species_id_.size();
    for (int j=0; j<num_states; j++) {
      for (int i=0; i<num_states; i++) {
        dense_id = num_states*j + i;

        //Dense rows and columns for local mdot
        if(j==num_states-2 || i==num_states-2) {
          //dense_to_sparse_[dense_id] = 1; //already in getjacobianpattern
        }
        //Dense rows and columns for local T
        if(j==num_states-1 || i==num_states-1) {
          //dense_to_sparse_[dense_id] = 1; //already in getjaocobian patter
        }
        //Dense rows for off-diagonal T?
        if(i==num_states-1) {
          //dense_to_sparse_offdiag_[dense_id] = 1;
        }

        //Diagonal
        if(i==j) {
          dense_to_sparse_[dense_id] = 1;
          dense_to_sparse_offdiag_[dense_id] = 1;
        }
        //Dense columns? for fuel species
                for (int k=0; k<num_fuel_species; k++) {
          //if(i==fuel_species_id_[k] || j==fuel_species_id_[k]) {
          if(j==fuel_species_id_[k]) {
            //dense_to_sparse_[dense_id] = 1;
            //dense_to_sparse_offdiag_[dense_id] = 1;
          }
        }
        //Dense columns? for oxidizer species
	for (int k=0; k<num_oxidizer_species; k++) {
          //if(i==oxidizer_species_id_[k] || j==oxidizer_species_id_[k]) {
          if(j==oxidizer_species_id_[k]) {
            //dense_to_sparse_[dense_id] = 1;
            //dense_to_sparse_offdiag_[dense_id] = 1;
          }
	}
      }
    }

    // Count non-zeros
    int width = 2*num_off_diagonals_ + 1;
    int i1, i2;
    for (int group=1; group<=width; group++) {
      for (int j=group-1; j<num_local_points*num_states; j+=width) {
        int jglobal = j + my_pe_*num_local_points*num_states;
        int jstate = jglobal % num_states;
        i1 = max(0, jglobal - jstate - num_states);
        i2 = min(jglobal + (num_states-1 - jstate) + num_states, num_points*num_states-1);
	for (int i=i1; i<=i2; i++) {
          int istate = i % num_states;
          dense_id = num_states*jstate + istate; //here j is column and i is row
          //Diagonal block
          if (i>= jglobal-jstate && i<=jglobal+num_states-1-jstate) {
            if (dense_to_sparse_[dense_id] == 1)
              num_nonzeros_loc_++;
          }
          //Off-diagonal block
          if (i<jglobal-jstate || i>jglobal+num_states-1-jstate) {
            if (dense_to_sparse_offdiag_[dense_id] == 1)
              num_nonzeros_loc_++;
          }
        }
      }
    }

    // SuperLU
    reactor_jacobian_dist_.assign(num_nonzeros_loc_, 0.0);
    col_id_.assign(num_nonzeros_loc_, 0);
    row_id.assign(num_nonzeros_loc_, 0);
    row_sum_.assign(num_local_points*num_states+1,0);

    // Get pattern "manually" for now
    int innz=0;
    // For compressed row storage j is the row and i is the column
    for (int j=0; j<num_local_points*num_states; j++) {
      int jglobal = j + my_pe_*num_local_points*num_states;
      int jstate = jglobal % num_states;
      i1 = max(0, jglobal - jstate - num_states);
      i2 = min(jglobal + (num_states-1 - jstate) + num_states, num_points*num_states-1);
      for (int i=i1; i<=i2; i++) {
        // Compressed column storage
        // row_id_[innz] = i;
	// column_id[innz] = j;

        // Compressed row storage
        int istate = i % num_states;
        dense_id = num_states*istate + jstate; //i is column and j is row

        //Diagonal block. Should use GetJacobianPattern?
        if (i>= jglobal-jstate && i<=jglobal+num_states-1-jstate) {
          if (dense_to_sparse_[dense_id] == 1) {
            row_id[innz] = j;
            col_id_[innz] = i;
            innz++;
          }
        }
        //Off-diagonal block.
        if (i<jglobal-jstate || i>jglobal+num_states-1-jstate) {
          if (dense_to_sparse_offdiag_[dense_id] == 1) {
            row_id[innz] = j;
            col_id_[innz] = i;
            innz++;
          }
        }
      } // for i=i1 to i2
    } //for j=0 to < num_local_states

    for(int j=0; j<num_nonzeros_loc_; ++j) {
        ++row_sum_[row_id[j]+1];
    } // for(int j=0; j<num_nonzeros; ++j)

    // At this point the column sum contains just the count per column, need
    // to loop through to change the count to a sum
    for(int j=0; j<num_local_points*num_states; ++j) {
      row_sum_[j+1] += row_sum_[j];
    }

  } // if (integrator_type_ == 2)

 // Set sparse matrix
  if(integrator_type_ == 2) {
    if(superlu_serial_) {
      sparse_matrix_ = new SparseMatrix(num_local_points*num_states, num_nonzeros_loc_);
      if(sparse_matrix_ == NULL) {
        printf("fail sparse_matrix_ \n");
        logger_->PrintF("# ERROR: failed to allocate new SparseMatrix object\n");
        exit(-1); // TODO: add recoverable failure
      }
    } else {
      sparse_matrix_dist_ = new SparseMatrix_dist(num_local_points*num_states, num_nonzeros_loc_,comm_);
      if(sparse_matrix_dist_ == NULL) {
        printf("fail sparse_matrix_dist_ \n");
        logger_->PrintF("# ERROR: failed to allocate new SparseMatrix object\n");
        exit(-1); // TODO: add recoverable failure
      }
    }

    if(store_jacobian_) {
      // TODO: add allocation check in case we run out of memory for large
      //       mechanisms and grids
      //saved_jacobian_.assign(num_nonzeros, 0.0);
      saved_jacobian_dist_.assign(num_nonzeros_loc_, 0.0);
    }
  } // if integrator_type == 2

  // Case 3 is an approximately factorized matrix
  // local sparse matrices at each grid point for the chemical Jacobian (SuperLU)
  // global tridiagonal matrix for the transport Jacobian (LAPACK)
  if(integrator_type_ == 3) {
    //SuperLU
    row_id_chem_.assign(num_nonzeros_zerod, 0);
    column_id_chem_.assign(num_nonzeros_zerod, 0);
    column_sum_chem_.assign(num_states+1,0);
    reactor_jacobian_chem_.assign(num_nonzeros_zerod, 0.0);

    // The Jacobian pattern is assumed to be in compressed column storage
    reactor_->GetJacobianPattern(&row_id_chem_[0],
                                 &column_id_chem_[0]);

    for(int j=0; j<num_nonzeros_zerod; ++j) {
      ++column_sum_chem_[column_id_chem_[j]+1];
    }

    for(int j=0; j<num_states; ++j) {
      column_sum_chem_[j+1] += column_sum_chem_[j];
    }


    diagonal_id_chem_.assign(num_states, -1); // use negative one as flag to
    // check if there are any missing
    // diagonal terms that would prevent
    // its use for the sparse factorization
    for(int j=0; j<num_nonzeros_zerod; ++j) {
      if(row_id_chem_[j] == column_id_chem_[j]) {
        diagonal_id_chem_[row_id_chem_[j]] = j;
      }
    }
    for(int j=0; j<num_states; ++j) {
      if(diagonal_id_chem_[j] == -1) {
        logger_->PrintF(
                        "# ERROR: reactor Jacobian does not have diagonal in row %d\n",
                        j);

        valid_jacobian_structure_ = false;
        exit(-1); // TODO: add recoverable failure
      }
   }

    sparse_matrix_chem_.assign(num_reactors, NULL);
    for(int j=0; j<num_reactors; ++j) {
      sparse_matrix_chem_[j] = new SparseMatrix(num_states, num_nonzeros_zerod);
      if(sparse_matrix_chem_[j] == NULL) {
	logger_->PrintF(
          "# ERROR: failed to allocate new SparseMatrix object\n"
          "#        for grid point %d (z = %.18g [m])\n",
          j, z_[j]);
	exit(-1); // TODO: add recoverable failure
      }
    }

    if(store_jacobian_) {
      saved_jacobian_chem_.assign(num_nonzeros_zerod*num_reactors, 0.0);
    }

    // LAPACK
    int storage = 4*1 + 1;
    banded_jacobian_.assign(num_local_points*storage * num_states, 0.0);

    // test over whole domain
    num_states_per_proc_ = (num_states + npes_ - 1)/npes_;
    if(num_states - (npes_-1)*num_states_per_proc_ < 1) {
      //num_states_per_proc_ -= 1;
      //if(num_states_per_proc_ < 1) {
        cerr << "Too few states per processor. Try using fewer processors.\n";
        exit(-1);
        //}
    }
    if(my_pe_ == npes_-1) {
      num_states_local_ = num_states - (npes_-1)*num_states_per_proc_;
    } else {
      num_states_local_ = num_states_per_proc_;
    }
    banded_jacobian2_.assign(num_points*storage * num_states_local_, 0.0);
    banded_jacobian_serial_.assign(num_points*4*num_states_local_, 0.0);
    pivots_serial_.assign(num_points*num_states_local_, 0);
  } // if integrator_type == 3

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
