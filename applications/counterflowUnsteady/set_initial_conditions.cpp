#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
#include <utilities/string_utilities.h>
#include <utilities/math_utilities.h>

#include "set_initial_conditions.h"



// Set initial composition
void SetInitialComposition(FlameParams &flame_params, double *y, double *time)
{
  const int num_points = flame_params.z_.size();
  const int num_local_points  = flame_params.num_local_points_;
  const int num_states  = flame_params.reactor_->GetNumStates();
  const int num_species = flame_params.reactor_->GetNumSpecies();

  const double pressure = flame_params.parser_->pressure();
  const double ref_temperature = flame_params.parser_->ref_temperature();
  const double ref_momentum = flame_params.parser_->ref_momentum();

  const int nover = flame_params.nover_;

  int my_pe = flame_params.my_pe_;
  MPI_Comm comm = flame_params.comm_;

  double species_left, species_right, species_stoich, species;
  double temperature_left, temperature_right, temperature,
    temperature_stoich, temperature_eq;
  double velocity_fuel, velocity_oxidizer;
  double thickness;
  double relative_volume, relative_volume_eq;
  double mixture_molecular_mass = 0.0;
  std::vector<double> species_molecular_mass;
  std::vector<double> eq_mass_fractions;
  double length = flame_params.length_;
  double strain_rate, x0, f, conductivity, cp, zeta, zmix, zst;
  double sumY;
  int transport_error;

  double beta = 2.0; //1.0 + flame_params.simulation_type_; // only axisymmetric

  std::string state_name= "MassFraction_OH";
  int OH_id = flame_params.reactor_->GetIdOfState(state_name.c_str());

  species_molecular_mass.assign(num_species, 0.0);
  flame_params.reactor_->GetSpeciesMolecularWeight(&species_molecular_mass[0]);

  eq_mass_fractions.assign(num_species, 0.0);

  temperature_stoich = flame_params.parser_->ignition_temperature()/ref_temperature;

  x0 = length*0.5; //put flame in the middle

  double centerWidth = flame_params.parser_->thickness();
  double slopeWidth = flame_params.parser_->thickness();
  double scale = 0.8*length / (centerWidth + 2*slopeWidth);
  if(scale < 1.0) {
    centerWidth *= scale;
    slopeWidth *= scale;
  }
  double dz = flame_params.z_[int(num_points/2)] - flame_params.z_[int(num_points/2)-1];
  int centerPointCount = int(0.5 + 0.5*centerWidth/dz);
  int slopePointCount = int(0.5 + slopeWidth/dz);
  int jm = int(num_points /2);
  int jl2 = jm - centerPointCount;
  int jl1 = jl2 - slopePointCount;
  int jr1 = jm + centerPointCount;
  int jr2 = jr1 + slopePointCount;

  if(flame_params.flame_type_ == 1) {
    x0 *= 2;
    jl1 *= 2;
    jl2 *= 2;
    jr1 *= 2;
    jr2 *= 2;
  }

  zst = flame_params.stoichiometric_mixture_fraction_;

  /*
  // ------------ BEGIN Constrained Equibrium calc  -----------//
  // Only on root
  if(my_pe == 0) {
    printf("# Using equilibrium composition to initialize Y and T\n");
    int k,l;
    std::vector<double> thermoCoeffs(num_species*16);
    std::vector<double> thermoCoeffs_CEQ(num_species*15);
    flame_params.mechanism_->getThermoCoeffs(&thermoCoeffs[0]);
    for(k=0; k<num_species; k++)
      for(l=0; l<15; l++) //transpose
        thermoCoeffs_CEQ[l*num_species+k] = thermoCoeffs[k*16 + l];

    int ne;
    ne = flame_params.mechanism_->getNumElements();
    std::vector<int> Ein_int(num_species*ne);
    std::vector<double> Ein_double(num_species*ne);
    flame_params.mechanism_->getSpeciesOxygenCount(&Ein_int[0]);
    flame_params.mechanism_->getSpeciesNitrogenCount(&Ein_int[2*num_species]);
    flame_params.mechanism_->getSpeciesHydrogenCount(&Ein_int[num_species]);
    if(ne>3) flame_params.mechanism_->getSpeciesCarbonCount(&Ein_int[3*num_species]);
    if(ne>4) flame_params.mechanism_->getSpeciesArgonCount(&Ein_int[4*num_species]);
    if(ne>5) flame_params.mechanism_->getSpeciesHeliumCount(&Ein_int[5*num_species]);

    for(k=0; k<num_species*ne; k++)
      Ein_double[k] = (double)Ein_int[k];

    ofstream file1;
    file1.open("CEQ-inputs");
    file1 << num_species << "\n";
    file1 << ne << "\n";
    file1 << pressure << "\n";
    if(flame_params.flame_type_ == 0) {
      file1 << 0.5*(flame_params.fuel_temperature_ +
                    flame_params.oxidizer_temperature_)*
        flame_params.ref_temperature_ << "\n";
    } else {
      file1 << flame_params.fuel_temperature_*flame_params.ref_temperature_ << "\n";
    }
    for(k=0; k<num_species; k++)
      file1 << 1.0/flame_params.inv_molecular_mass_[k] << " ";
    file1 << "\n";
    for(k=0; k<num_species*15; k++)
      file1 << thermoCoeffs_CEQ[k] << " ";
    file1 << "\n";
    for(k=0; k<num_species*ne; k++)
      file1 << Ein_double[k] << " ";
    file1 << "\n";
    for(k=0; k<num_species; k++) {
      if(flame_params.flame_type_ == 0) {
        file1 << zst*flame_params.fuel_mass_fractions_[k] + (1.0-zst)*flame_params.oxidizer_mass_fractions_[k] << " ";
      } else {
        file1 << flame_params.inlet_mass_fractions_[k] << " ";
      }
    }
    file1 << "\n";
    file1.close();

    // Call CEQ
    system("/usr/apps/advcomb/bin/eqHPfromFile.x");

    // Read equilibrium state
    std::ifstream infile("CEQ.dat");
    std::string line;
    k=0;
    double val;
    while (std::getline(infile, line))
    {
      std::istringstream iss(line);
      iss >> val;
      if(k<num_species)
        eq_mass_fractions[k] = val;
      if(k==num_species)
        temperature_eq = val/flame_params.ref_temperature_;

      k++;
    }
    printf("# Eq. temp: %5.3e\n",temperature_eq);
  }

  // Broadcast eq. T & Y to all procs
  MPI_Bcast(&temperature_eq, 1, MPI_DOUBLE, 0, comm);
  MPI_Bcast(&eq_mass_fractions[0], num_species, MPI_DOUBLE, 0, comm);
  // ------------ END Constrained Equibrium calc  -----------//
  */

  /**/
  // "Simple" initialization sometimes works better than EQ
  sumY = 0.0;
  for(int k=0; k<num_species; ++k){
    if(flame_params.flame_type_ == 0) {
      eq_mass_fractions[k] = zst*flame_params.fuel_mass_fractions_[k] +
        (1.0-zst)*flame_params.oxidizer_mass_fractions_[k];
      } else {
      eq_mass_fractions[k] = flame_params.inlet_mass_fractions_[k];
      }

    if(k==OH_id)
      eq_mass_fractions[k] += 0.003; // add OH to help ignition
    sumY += eq_mass_fractions[k];
  }
  for(int k=0; k<num_species; ++k)
    eq_mass_fractions[k] /= sumY;

  temperature_eq = temperature_stoich;
  /**/

  // Renormalize and compute mixture weight
  mixture_molecular_mass = 0.0;
  for(int k=0; k<num_species; ++k) {
    mixture_molecular_mass += eq_mass_fractions[k]/species_molecular_mass[k];
  }
  mixture_molecular_mass = 1.0/mixture_molecular_mass;

  // Relative volume
  relative_volume_eq = flame_params.reactor_->GetGasConstant()*
    ref_temperature*temperature_eq/(pressure*mixture_molecular_mass);


  // Compute inlet velocities to approximate strain rate
  if(flame_params.flame_type_ == 0) {
    velocity_fuel = flame_params.mass_flux_fuel_*flame_params.fuel_relative_volume_;
    velocity_oxidizer = flame_params.mass_flux_oxidizer_*flame_params.oxidizer_relative_volume_;
  } else if (flame_params.flame_type_ == 1) {
    velocity_fuel = flame_params.mass_flux_fuel_*flame_params.inlet_relative_volume_;
    velocity_oxidizer = 0.0;
  } else if (flame_params.flame_type_ == 2) {
    velocity_fuel = flame_params.mass_flux_fuel_*flame_params.inlet_relative_volume_;
    velocity_oxidizer = flame_params.mass_flux_oxidizer_*relative_volume_eq;
  }

  if(flame_params.parser_->finite_separation()) {
    strain_rate = (fabs(velocity_fuel) +
                   fabs(velocity_oxidizer))/length;
  } else {
    strain_rate = flame_params.strain_rate_;
  }

  // For R-to-P flame, overwrite oxidizer quantities with eq. quantities
  if(flame_params.flame_type_ == 2) {
    flame_params.oxidizer_temperature_ = temperature_eq;
    flame_params.oxidizer_relative_volume_ = relative_volume_eq;
    for(int k=0; k<num_species; ++k) {
      flame_params.oxidizer_mass_fractions_[k] = eq_mass_fractions[k];
    }
  }

  // Set right/left T
  if(flame_params.flame_type_ == 0 or flame_params.flame_type_ == 2) {
    temperature_right = flame_params.oxidizer_temperature_;
    temperature_left = flame_params.fuel_temperature_;
  } else if (flame_params.flame_type_ == 1) {
    temperature_left = flame_params.fuel_temperature_;
    temperature_right= flame_params.fuel_temperature_;
  }

  for(int j=0; j<num_local_points; ++j) {
    int jglobal = j + my_pe*num_local_points;
    int jext = j + nover;

    double ramp = 0.0;
    if(jglobal <= jl1) {
      temperature = temperature_left;
      for(int k=0; k<num_species; ++k) {
        if(flame_params.flame_type_ == 0) {
          y[j*num_states+k] = flame_params.fuel_mass_fractions_[k];
        } else {
          y[j*num_states+k] = flame_params.inlet_mass_fractions_[k];
        }
      }
    } else if (jglobal > jl1 and jglobal <= jl2) {
      ramp = ((double)jglobal-(double)jl1)/((double)jl2-(double)jl1);
      temperature = temperature_left + (temperature_eq-temperature_left)*ramp;
      for(int k=0; k<num_species; ++k) {
        if(flame_params.flame_type_ == 0) {
          y[j*num_states+k] = flame_params.fuel_mass_fractions_[k] +
            (eq_mass_fractions[k]-flame_params.fuel_mass_fractions_[k])*ramp;
        } else {
          y[j*num_states+k] = flame_params.inlet_mass_fractions_[k] +
            (eq_mass_fractions[k]-flame_params.inlet_mass_fractions_[k])*ramp;
        }
      }
    } else if (jglobal > jl2 and jglobal <= jr1) {
      temperature = temperature_eq;
      for(int k=0; k<num_species; ++k) {
        y[j*num_states+k] = eq_mass_fractions[k];}
    } else if (jglobal > jr1 and jglobal <= jr2) {
      ramp = ((double)jglobal-(double)jr1)/((double)jr2-(double)jr1);
      temperature = temperature_eq - (temperature_eq-temperature_right)*ramp;
      for(int k=0; k<num_species; ++k) {
        if(flame_params.flame_type_ == 0) {
          y[j*num_states+k] = eq_mass_fractions[k] -
            (eq_mass_fractions[k]-flame_params.oxidizer_mass_fractions_[k])*ramp;
        } else {
          y[j*num_states+k] = eq_mass_fractions[k] -
            (eq_mass_fractions[k]-flame_params.inlet_mass_fractions_[k])*ramp;
        }
      }
    } else {
      temperature = temperature_right;
      for(int k=0; k<num_species; ++k) {
        if(flame_params.flame_type_ == 0) {
          y[j*num_states+k] = flame_params.oxidizer_mass_fractions_[k];
        } else {
          y[j*num_states+k] = flame_params.inlet_mass_fractions_[k];
        }
      }
    }

    // Renormalize and compute mixture weight
    sumY = 0.0;
    for(int k=0; k<num_species; ++k)
      sumY += y[j*num_states+k];

    mixture_molecular_mass = 0.0;
    for(int k=0; k<num_species; ++k) {
      y[j*num_states+k] /= sumY;
      mixture_molecular_mass += y[j*num_states+k]/species_molecular_mass[k];
    }
    mixture_molecular_mass = 1.0/mixture_molecular_mass;

    // Relative volume
    relative_volume = flame_params.reactor_->GetGasConstant()*
      ref_temperature*temperature/(pressure*mixture_molecular_mass);

    y[j*num_states + num_species] = relative_volume;

    // Mass flux
    if(flame_params.parser_->finite_separation()) {
      flame_params.mass_flux_[j] = velocity_fuel + flame_params.z_[jglobal]/length*
        (-fabs(velocity_oxidizer)-velocity_fuel);
      flame_params.mass_flux_[j] /= relative_volume;
    }

    // Temperature
    y[j*num_states + num_species+1] = temperature;

    // Momentum
    if(flame_params.parser_->finite_separation()) {
      if(flame_params.z_[jglobal] <= x0)
        y[j*num_states + num_species+2] = flame_params.z_[jglobal]/x0*strain_rate/ref_momentum;
      else
        y[j*num_states + num_species+2] =
          (strain_rate-(flame_params.z_[jglobal]-x0)/(length-x0)*strain_rate)/ref_momentum;
    } else {
      y[j*num_states + num_species+2] = flame_params.strain_rate_/beta*
        sqrt(relative_volume/flame_params.oxidizer_relative_volume_)/
        ref_momentum;
    }

    if(flame_params.parser_->finite_separation()) {
      // Pstrain/lambda -- initialize to a^2/beta^2
      y[j*num_states + num_species + 3] = -strain_rate*strain_rate/beta/beta/ref_momentum;
    }

  } // j spatial loop

  // gather y to have all grid points
  if(!flame_params.parser_->finite_separation()) {
    std::vector<double> mass_flux;
    mass_flux.assign(num_points,0.0);

    long int dsize = num_local_points*num_states;
    double *ybuf;
    if(my_pe == 0)
      ybuf = (double *)malloc(flame_params.npes_*dsize*sizeof(double));

    int nodeDest, nodeFrom;
    nodeDest = 0; // gather to root
    nodeFrom = 0; // scatter from root

    MPI_Gather(&y[0],
               dsize,
               PVEC_REAL_MPI_TYPE,
               ybuf,
               dsize,
               PVEC_REAL_MPI_TYPE,
               nodeDest,
               flame_params.comm_);

    if(my_pe == 0){
      // Mass flux
      int stag;
      if(flame_params.flame_type_ == 0 or flame_params.flame_type_ == 2) {
        stag = num_points/4;
      } else if (flame_params.flame_type_ == 1) {
        stag = num_points-1;
      }

      mass_flux[stag] = 0.0;
      for(int j=stag+1; j<num_points; j++) {
        mass_flux[j] = mass_flux[j-1] -
          ybuf[j*num_states + num_species+2]/ybuf[j*num_states + num_species]*
          beta*flame_params.dz_[j]*ref_momentum;
      }

      for(int j=stag-1; j>-1; j--) {
        mass_flux[j] = mass_flux[j+1] +
          ybuf[j*num_states + num_species+2]/ybuf[j*num_states + num_species]*
          beta*flame_params.dz_[j+1]*ref_momentum;
      }
    }

    dsize = num_local_points;
    MPI_Scatter(&mass_flux[0],
                dsize,
                PVEC_REAL_MPI_TYPE,
                &flame_params.mass_flux_[0],
                dsize,
                PVEC_REAL_MPI_TYPE,
                nodeFrom,
                flame_params.comm_);
  }

  // Local mass flux vector
  for(int j=0; j<num_local_points; j++) {
    int jext = j + nover;
    flame_params.mass_flux_ext_[jext] =  flame_params.mass_flux_[j];
  }

  // Restart using binary file if provided
  if(flame_params.parser_->restart_file() != std::string("/dev/null")) {
    const char* filename1;//[32];
    char filename2[32];
    int num_points_file, num_vars_file;

    int npes;
    int disp;
    std::vector<double> buffer;
    buffer.assign(num_local_points, 0.0);
    double time_file = 0.0;

    MPI_File restart_file;
    npes = flame_params.npes_;

    filename1 = flame_params.parser_->restart_file().c_str();
    sprintf(filename2,"%s",filename1);

    MPI_File_open(MPI_COMM_WORLD, filename2, MPI_MODE_RDONLY, MPI_INFO_NULL, &restart_file);

    MPI_File_read_all(restart_file, &num_points_file, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_all(restart_file, &num_vars_file, 1, MPI_INT, MPI_STATUS_IGNORE);
    if(num_vars_file != num_states) {
      cerr << "WARNING: restart file and mechanism have different number of species. Species not found will be initialized at 0.\n";
    }

    MPI_File_read_all(restart_file, &time_file, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
    *time = time_file;

    string file_state_names[num_vars_file];
    for(int j=0; j<num_vars_file; ++j) {
      char buf[64];
      MPI_File_read_all(restart_file, &buf, 64, MPI_CHAR, MPI_STATUS_IGNORE);
      file_state_names[j] = string(buf);
    }

    // Initialize y to 0
    for (int k=0; k<num_local_points*num_states; ++k)
      y[k] = 0.0;

    for(int j=0; j<num_states; ++j) {
      string state_name = flame_params.reactor_->GetNameOfStateId(j);
      for(int i=0; i<num_vars_file; ++i) {
        //if(state_name == file_state_names[i]) {
        if(strcasecmp(state_name.c_str(),file_state_names[i].c_str()) == 0 ) {
          // Read restart_file
          disp = 2*sizeof(int) + sizeof(double) + num_vars_file*sizeof(char)*64
            + i*2*sizeof(double) // Left & right BCs from previous variables
            + sizeof(double) // Left BC from current variable
            + (my_pe + i*npes)*num_local_points*sizeof(double);
          MPI_File_set_view(restart_file, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
          MPI_File_read(restart_file, &buffer[0], num_local_points, MPI_DOUBLE, MPI_STATUS_IGNORE);
          // Set y values
          for (int k=0; k<num_local_points; ++k) {
            y[k*num_states + j] = buffer[k];
            if(j==num_species+1){
              y[k*num_states + j] /= flame_params.ref_temperature_;
            }
            if(j==num_species+2){
              y[k*num_states + j] /= flame_params.ref_momentum_;
            }
            if(j==num_species+3){
              y[k*num_states + j] /= flame_params.ref_momentum_;
            }
          }
          // Exit i loop
          break;
        }
        if (i==num_vars_file-1 && state_name != file_state_names[i]) {
          // Variable not found
          cerr << "WARNING: " << state_name << " not found in restart file.\n";
        }
      } // for i<num_vars_file
    } // for j<num_states

    // Read mass flux
    disp = 2*sizeof(int) + sizeof(double) + num_vars_file*sizeof(char)*64
      + num_vars_file*2*sizeof(double) // Left & right BCs from previous variables
      + sizeof(double) // Left BC from current variable
      + (my_pe + num_vars_file*npes)*num_local_points*sizeof(double);
    MPI_File_set_view(restart_file, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    MPI_File_read(restart_file, &buffer[0], num_local_points, MPI_DOUBLE, MPI_STATUS_IGNORE);

    for (int k=0; k<num_local_points; ++k) {
      flame_params.mass_flux_[k] = buffer[k];
      flame_params.mass_flux_ext_[k+2] = flame_params.mass_flux_[k];
    }

    MPI_File_close(&restart_file);

  } // if restart_file

}
