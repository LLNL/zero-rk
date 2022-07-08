#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
#include <utilities/string_utilities.h>
#include <utilities/math_utilities.h>
#include <utilities/file_utilities.h>

#include "set_initial_conditions.h"


// Set the mass fractions and temperature to the fuel composition
// if the grid point position is less than or equal to the
// stoichiometric mixture fraction and to the oxidizer composition
// otherwise.
void SetInitialCompositionAndWallTemp(FlameParams &flame_params, double *y, double *time)
{
  //  const int num_points  = flame_params.z_.size();
  const int num_local_points  = flame_params.num_local_points_;
  const int num_states  = flame_params.reactor_->GetNumStates();
  const int num_species = flame_params.fuel_mass_fractions_.size();

  const double pressure = flame_params.parser_->pressure();
  const double ref_temperature = flame_params.parser_->ref_temperature();

  int my_pe = flame_params.my_pe_;

  double species_left, species_right, delta_species, species;
  double temperature_left, temperature_right, delta_temperature_left,
    delta_temperature_right, temperature, temperature_stoich;
  double thickness;
  double mixture_molecular_mass = 0.0;
  std::vector<double> species_molecular_mass;

  // Add OH to facilitate ignition
  std::string state_name= "MassFraction_OH";
  int OH_id = flame_params.reactor_->GetIdOfState(state_name.c_str());

  species_molecular_mass.assign(num_species, 0.0);
  flame_params.reactor_->GetSpeciesMolecularWeight(&species_molecular_mass[0]);

  temperature_right = flame_params.fuel_temperature_;
  temperature_left = flame_params.oxidizer_temperature_;
  temperature_stoich = flame_params.parser_->ignition_temperature()/ref_temperature;
  delta_temperature_left = temperature_stoich - temperature_left;
  delta_temperature_right = temperature_stoich - temperature_right;

  thickness = flame_params.parser_->thickness();

  for(int j=0; j<num_local_points; ++j) {
    int jglobal = j + my_pe*num_local_points;

    if(flame_params.fix_temperature_) {
      temperature = flame_params.fixed_temperature_[jglobal];
    } else {
      if(flame_params.z_[jglobal] <= flame_params.stoichiometric_mixture_fraction_) {
        temperature = temperature_left + delta_temperature_left*flame_params.z_[jglobal]/
          flame_params.stoichiometric_mixture_fraction_;
      } else {
        temperature = temperature_stoich - delta_temperature_right*
          (flame_params.z_[jglobal]-flame_params.stoichiometric_mixture_fraction_)/
          (1.0-flame_params.stoichiometric_mixture_fraction_);
      }
    }

    for(int k=0; k<num_species; ++k) {
      species_right = flame_params.fuel_mass_fractions_[k];
      species_left= flame_params.oxidizer_mass_fractions_[k];
      delta_species = species_right - species_left;

      species = (species_left + 0.5*delta_species) +
        0.5*delta_species*tanh(0.5*(flame_params.z_[jglobal]-flame_params.stoichiometric_mixture_fraction_)/thickness);

      y[j*num_states+k] = species;
    }

    // ADD OH
    if(abs(flame_params.z_[jglobal] -flame_params.stoichiometric_mixture_fraction_) <= 2.0*thickness) {
      y[j*num_states + OH_id] += 0.003;
    }

    // Renormalize and compute mixture weight
    double sumY = 0.0;
    for(int k=0; k<num_species; ++k)
      sumY += y[j*num_states+k];

    mixture_molecular_mass = 0.0;
    for(int k=0; k<num_species; ++k) {
      y[j*num_states+k] /= sumY;
      mixture_molecular_mass += y[j*num_states+k]/species_molecular_mass[k];
    }
    mixture_molecular_mass = 1.0/mixture_molecular_mass;

    y[j*num_states + num_species]   = flame_params.reactor_->GetGasConstant()*ref_temperature*
      temperature/(pressure*mixture_molecular_mass);

    y[j*num_states + num_species + 1] = temperature;

  } // j spatial loop


  // Restart using binary file if provided
  if(flame_params.parser_->restart_file() != std::string(zerork::utilities::null_filename)) {
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

    std::vector<string> file_state_names(num_vars_file);
    for(int j=0; j<num_vars_file; ++j) {
      char buf[64];
      MPI_File_read_all(restart_file, &buf, 64, MPI_CHAR, MPI_STATUS_IGNORE);
      file_state_names[j] = string(buf);
    }

    // Initialize y to 0
    for (int k=0; k<num_local_points*num_states; ++k) {
      y[k] = 0.0;
    }

    /**/
    for(int j=0; j<num_states; ++j) {
      string state_name = zerork::utilities::GetLowerCase(flame_params.reactor_->GetNameOfStateId(j));
      for(int i=0; i<num_vars_file; ++i) {
        string file_state_name = zerork::utilities::GetLowerCase(file_state_names[i]);
        if(state_name == file_state_name) {
          // Skip over BC data

          // Read interior data
          disp = 2*sizeof(int) + sizeof(double) + num_vars_file*sizeof(char)*64
            + i*2*sizeof(double) // Left & right BCs from previous variables
            + sizeof(double) // Left BC from current variable
            + (my_pe + i*npes)*num_local_points*sizeof(double);
          MPI_File_set_view(restart_file, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
          MPI_File_read(restart_file, &buffer[0], num_local_points, MPI_DOUBLE, MPI_STATUS_IGNORE);
          // Set y values
          for (int k=0; k<num_local_points; ++k) {
            y[k*num_states + j] = buffer[k];
            if(j==num_states-1){
              y[k*num_states + j] /= flame_params.ref_temperature_;
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
    /**/

    /*
    for(int j=0; j<num_states; ++j) {
      int i=j;
      // Read interior data
      disp = 2*sizeof(int) + sizeof(double) + num_vars_file*sizeof(char)*64
        + i*2*sizeof(double) // Left & right BCs from previous variables
        + sizeof(double) // Left BC from current variable
        + (my_pe + i*npes)*num_local_points*sizeof(double);
      MPI_File_set_view(restart_file, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
      MPI_File_read(restart_file, &buffer[0], num_local_points, MPI_DOUBLE, MPI_STATUS_IGNORE);
      // Set y values
      for (int k=0; k<num_local_points; ++k) {
        y[k*num_states + j] = buffer[k];
        if(j==num_states-1){
          y[k*num_states + j] /= flame_params.ref_temperature_;
        }
      }
    } // for j<num_states
    */

    MPI_File_close(&restart_file);

  } // if restart_file

}
