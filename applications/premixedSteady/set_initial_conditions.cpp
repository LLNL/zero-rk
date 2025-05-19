#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
#include <utilities/string_utilities.h>
#include <utilities/math_utilities.h>
#include <utilities/file_utilities.h>

#include "set_initial_conditions.h"

#include "equil/zerork_equilibrium_solver.h"

// Set all grid points to the inlet conditions including temperature
void SetConstantInlet(FlameParams &flame_params, double *y)
{
  //  const int num_points  = flame_params.z_.size();
  const int num_local_points  = flame_params.num_local_points_;
  const int num_states  = flame_params.reactor_->GetNumStates();
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

  for(int j=0; j<num_local_points; ++j) {

    for(int k=0; k<num_species; ++k) {
      y[j*num_states+k] = flame_params.inlet_mass_fractions_[k];
    }
    y[j*num_states + num_species]   = relative_volume;
    y[j*num_states + num_species+1] = inlet_temperature;

  }

}

// Set initial composition and temperature
// If it's a flame speed calculation, species are set to inlet values throughout the domain.
// Temperature is set to inlet value in the left half and to ignition_temperature on the right.
// If it's a flame-in-tube problem, species are set to inlet values on the left and to
// initial values on the right. The temperature is set to the wall temperature.
// Those values are overwritten if a restart file is specified
void SetInitialCompositionAndWallTemp(FlameParams &flame_params, double *y, double *time)
{
  const int num_local_points  = flame_params.num_local_points_;
  const int num_states  = flame_params.reactor_->GetNumStates();
  const int num_species = flame_params.inlet_mass_fractions_.size();

  const double pressure = flame_params.parser_->pressure();
  const double ref_temperature = flame_params.parser_->ref_temperature();

  int my_pe = flame_params.my_pe_;

  double relative_volume, relative_volume_left, relative_volume_right, delta_relative_volume;
  double species_left, species_right, delta_species, species;
  double inlet_molecular_mass = 0.0;
  double initial_molecular_mass = 0.0;
  double temp_left, temp_right, delta_temp;
  std::vector<double> species_molecular_mass;

  species_molecular_mass.assign(num_species, 0.0);
  flame_params.reactor_->GetSpeciesMolecularWeight(&species_molecular_mass[0]);
  for(int k=0; k<num_species; ++k) {
    initial_molecular_mass +=
      flame_params.initial_mass_fractions_[k]/species_molecular_mass[k];
    inlet_molecular_mass +=
      flame_params.inlet_mass_fractions_[k]/species_molecular_mass[k];
  }
  initial_molecular_mass = 1.0/initial_molecular_mass;
  inlet_molecular_mass   = 1.0/inlet_molecular_mass;

  // For flame speed calculations
  if(flame_params.simulation_type_ == 1) {
#ifdef ZERORK_HAVE_EQUILIBRIUM_SOLVER
    if(flame_params.parser_->use_equilibrium_for_init()) {
      // ------------ BEGIN Constrained Equibrium calc  -----------//
      double equi_temperature;
      std::vector<double> equi_mass_fractions(num_species);

      // Only on root
      if(my_pe == 0) {
        printf("# Using equilibrium composition to initialize Y and T\n");

        equi_temperature = flame_params.inlet_temperature_*ref_temperature;
        for(int k=0; k<num_species; k++) {
          equi_mass_fractions[k] = flame_params.inlet_mass_fractions_[k];
        }
        zerork::equilibrium_solver zeq(*flame_params.mechanism_);
        zeq.equilibrate(pressure, equi_temperature,
                        &equi_mass_fractions[0]);
        equi_temperature /= ref_temperature;
        printf("# Eq. temp: %5.3e\n",equi_temperature);
      }

      // Broadcast eq. T & Y to all procs
#ifdef ZERORK_MPI
      MPI_Bcast(&equi_temperature, 1, MPI_DOUBLE, 0, flame_params.comm_);
      MPI_Bcast(&equi_mass_fractions[0], num_species, MPI_DOUBLE, 0, flame_params.comm_);
#endif
      // ------------ END Constrained Equibrium calc  -----------//
      /**/

      // Relative positions where to begin/end linear ramp
      double ramp_start = flame_params.parser_->initial_inlet_extent();
      double ramp_end = ramp_start + flame_params.parser_->thickness();

      for(int j=0; j<num_local_points; j++) {
        int jglobal = j + my_pe*num_local_points;

        if (flame_params.z_[jglobal] <= ramp_start) {
          for(int k=0; k<num_species; k++) {
            y[j*num_states + k] = flame_params.inlet_mass_fractions_[k];
          }
          y[j*num_states + num_species + 1] = flame_params.inlet_temperature_;
        } else if (flame_params.z_[jglobal] > ramp_start && flame_params.z_[jglobal] <= ramp_end) {
          double frac = (flame_params.z_[jglobal]-ramp_start)/(ramp_end-ramp_start);
          for(int k=0; k<num_species; k++) {
            y[j*num_states + k] = frac*equi_mass_fractions[k] +
              (1.0-frac)*flame_params.inlet_mass_fractions_[k];
          }
          y[j*num_states + num_species + 1] = frac*equi_temperature +
            (1.0-frac)*flame_params.inlet_temperature_;
        } else {
          for(int k=0; k<num_species; k++) {
            y[j*num_states + k] = equi_mass_fractions[k];
          }
          y[j*num_states + num_species + 1] = equi_temperature;
        }

        double molecular_mass = 0.0;
        for(int k=0; k<num_species; k++) {
          molecular_mass += y[j*num_states + k]/species_molecular_mass[k];
        }
        molecular_mass = 1.0/molecular_mass;

        y[j*num_states + num_species] = flame_params.reactor_->GetGasConstant()*ref_temperature*
            y[j*num_states + num_species+1]/(pressure*molecular_mass);
      }
    } else
#else
    if(flame_params.parser_->use_equilibrium_for_init()) {
      if(my_pe == 0) {
        printf("No equilibrium interface enabled.  Please re-configure/re-build Zero-RK to enable CEQ or Cantera equilibrium interface.\n");
        printf("Falling back to \"simple\" initialization\n");
      }
    }
#endif
    {
      for(int j=0; j<num_local_points; ++j) {
        int jglobal = j + my_pe*num_local_points;
        temp_left = flame_params.parser_->inlet_temperature()/
          flame_params.parser_->ref_temperature();
        temp_right = flame_params.parser_->ignition_temperature()/
          flame_params.parser_->ref_temperature();
        delta_temp = temp_right - temp_left;

        y[j*num_states + num_species+1] = (temp_left+0.5*delta_temp) +
          0.5*delta_temp*tanh(0.5*(flame_params.z_[jglobal]-flame_params.parser_->initial_inlet_extent())/flame_params.parser_->thickness());

        y[j*num_states + num_species] = flame_params.reactor_->GetGasConstant()*ref_temperature*
          y[j*num_states + num_species+1]/(pressure*inlet_molecular_mass);

        for(int k=0; k<num_species; ++k) {
          y[j*num_states+k] = flame_params.inlet_mass_fractions_[k];
        }
      }
    }
    // For flame in tube calculations
  } else {

    for(int j=0; j<num_local_points; ++j) {
      int jglobal = j + my_pe*num_local_points;

      relative_volume_left = flame_params.reactor_->GetGasConstant()*ref_temperature*
	flame_params.wall_temperature_[jglobal]/(pressure*inlet_molecular_mass);

      relative_volume_right = flame_params.reactor_->GetGasConstant()*ref_temperature*
	flame_params.wall_temperature_[jglobal]/(pressure*initial_molecular_mass);

      delta_relative_volume = relative_volume_right - relative_volume_left;

      relative_volume = (relative_volume_left+0.5*delta_relative_volume) +
	0.5*delta_relative_volume*tanh(0.5*(flame_params.z_[jglobal]-flame_params.parser_->initial_inlet_extent())/flame_params.parser_->thickness());

      for(int k=0; k<num_species; ++k) {
	species_left = flame_params.inlet_mass_fractions_[k];
	species_right= flame_params.initial_mass_fractions_[k];
	delta_species = species_right - species_left;

	species= (species_left + 0.5*delta_species) +
	  0.5*delta_species*tanh(0.5*(flame_params.z_[jglobal]-flame_params.parser_->initial_inlet_extent())/flame_params.parser_->thickness());

	y[j*num_states+k] = species;
      }

      y[j*num_states + num_species]   = relative_volume;
      y[j*num_states + num_species+1] = flame_params.wall_temperature_[jglobal];

    } // j spatial loop

  } // if-else simulation_type == 1


  // If initial_profile_file is provided, overwrite and
  // interpolate species and temperature profiles from FlameMaster solution
  if(flame_params.parser_->initial_profile_file() != std::string(zerork::utilities::null_filename)) {
    std::vector<double> file_position, file_temperature, file_density, file_mflux, file_fields;
    std::ifstream initial_profile_file;
    std::string line;
    std::string delimiters = std::string(zerork::utilities::WHITESPACE) + ",";
    std::vector<std::string> fields;

    initial_profile_file.open(flame_params.parser_->initial_profile_file().c_str());

    if(initial_profile_file) {
      while(zerork::utilities::GetAnyLine(initial_profile_file,&line)) {
	// remove everything before the comment character
	//      zerork::utilities::SplitBeforeDelimiter(line,
	//				       "!",
	//				       true,
	//				       &line);
	zerork::utilities::SplitStringToVector(line,
						delimiters,
						&fields);
	// only record valid data lines
	//	if(fields.size() == 2+num_species) {
	if(zerork::utilities::StringIsDouble(fields[0]) &&
	   zerork::utilities::StringIsDouble(fields[1])) {

	  file_position.push_back(atof(fields[0].c_str()));
	  file_temperature.push_back(atof(fields[2].c_str()));
	  //	  file_density.push_back(atof(fields[num_species+num_species+3].c_str()));
	  file_mflux.push_back(atof(fields[1].c_str()));
	  for(int k=0; k<num_species; ++k) {
	    file_fields.push_back(atof(fields[3+k].c_str()));
	  }
	}// if stringisdouble
	//	}
      } // while(zerork::utilities::GetAnyLine(wall_temperature_file,&line))

    } else {
      printf("# ERROR: could not open initial_profile_file = %s\n",
	     flame_params.parser_->initial_profile_file().c_str());
      exit(-1);
    }//if file

    if(file_position.size() == 0) {
      printf("# ERROR: no position and temperature data found in\n"
	     "#        initial_profile_file = %s\n",
	     flame_params.parser_->initial_profile_file().c_str());
      exit(-1);
    } // if file_position.size

    zerork::utilities::SortVectors(&file_position, &file_temperature);

    zerork::utilities::InterpolationTable table(file_position.size(),
						 &file_position[0],
						 &file_temperature[0],
						 zerork::utilities::LINEAR,
						 false); // use_interpolation
    zerork::utilities::InterpolationTable table2(file_position.size(),
						  &file_position[0],
						  //&file_density[0],
						  &file_mflux[0],
						  zerork::utilities::LINEAR,
						  false); // use_interpolation
    for(int j=0; j<num_local_points; ++j) {
      int jglobal = j + my_pe*num_local_points;
      //      y[j*num_states + num_species]   = 1.0/table2.Interpolate(flame_params.z_[jglobal]);
      y[j*num_states + num_species]   = table2.Interpolate(flame_params.z_[jglobal]);
      y[j*num_states + num_species+1] = table.Interpolate(flame_params.z_[jglobal]);
      y[j*num_states + num_species+1] /= flame_params.parser_->ref_temperature();
    }

    for(int k=0; k<num_species; ++k) {
      std::vector<double> file_species;
      file_species.assign(file_position.size(),0.0);

      for (int j=0; j<(int)file_position.size(); ++j) {
	file_species[j] = file_fields[k + j*num_species];
      }
      zerork::utilities::InterpolationTable table(file_position.size(),
						   &file_position[0],
						   &file_species[0],
						   zerork::utilities::LINEAR,
						   false); // use_interpolation

      for(int j=0; j<num_local_points; ++j) {
	int jglobal = j + my_pe*num_local_points;
	y[j*num_states + k] = table.Interpolate(flame_params.z_[jglobal]);
      } // for j points

    } // for k species

  } // if initial_profile_file isn't zerork::utilities::null_filename

  // If binary restart_file is provided, overwrite with that
  if(flame_params.parser_->restart_file() != std::string(zerork::utilities::null_filename)) {
    const char* filename1;//[32];
    char filename2[32];
    int num_points_file, num_vars_file;

    int npes;
    int disp;
    std::vector<double> buffer;
    buffer.assign(num_local_points, 0.0);
    double time_file = 0.0;

#ifdef ZERORK_MPI
    MPI_File restart_file;
#endif
    npes = flame_params.npes_;

    filename1 = flame_params.parser_->restart_file().c_str();
    sprintf(filename2,"%s",filename1);
#ifdef ZERORK_MPI
    MPI_File_open(flame_params.comm_, filename2, MPI_MODE_RDONLY, MPI_INFO_NULL, &restart_file);
#else
    ifstream restart_file(filename2, ios::in | ios::binary);
#endif

#ifdef ZERORK_MPI
    MPI_File_read_all(restart_file, &num_points_file, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_all(restart_file, &num_vars_file, 1, MPI_INT, MPI_STATUS_IGNORE);
#else
    restart_file.read((char*)&num_points_file, sizeof(int));
    restart_file.read((char*)&num_vars_file, sizeof(int));
#endif
    if(num_vars_file != num_states) {
      printf("WARNING: restart file and mechanism have different number of species. Species not found will be initialized at 0.\n");
    }

#ifdef ZERORK_MPI
    MPI_File_read_all(restart_file, &time_file, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
#else
    restart_file.read((char*)&time_file, sizeof(double));
#endif
    *time = time_file;

    std::vector<string> file_state_names(num_vars_file);
    for(int j=0; j<num_vars_file; ++j) {
      char buf[64];
#ifdef ZERORK_MPI
      MPI_File_read_all(restart_file, &buf, 64, MPI_CHAR, MPI_STATUS_IGNORE);
#else
      restart_file.read((char*)&buf, 64*sizeof(char));
#endif
      file_state_names[j] = string(buf);
    }

    // Initialize y to 0
    for (int k=0; k<num_local_points*num_states; ++k) {
      y[k] = 0.0;
    }
    double bc;

    // Read data for each variable
    for(int j=0; j<num_states; ++j) {
      string state_name = zerork::utilities::GetLowerCase(flame_params.reactor_->GetNameOfStateId(j));
      for(int i=0; i<num_vars_file; ++i) {
	string file_state_name = zerork::utilities::GetLowerCase(file_state_names[i]);
        if(state_name == file_state_name) {
          // Read restart_file
#ifdef ZERORK_MPI
          disp = 2*sizeof(int) + sizeof(double) + num_vars_file*sizeof(char)*64
            + i*sizeof(double) //left BC from prev vars
            + sizeof(double) //left BC from cur var
            + (my_pe + i*npes)*num_local_points*sizeof(double);
          MPI_File_set_view(restart_file, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
          MPI_File_read(restart_file, &buffer[0], num_local_points, MPI_DOUBLE, MPI_STATUS_IGNORE);
#else
          restart_file.read((char*)&bc, sizeof(double)); //read/skip BC
          restart_file.read((char*)&buffer[0], num_local_points*sizeof(double));
#endif
          // Set y values
          for (int k=0; k<num_local_points; ++k) {
            y[k*num_states + j] = buffer[k];
            if(j==num_states-2){
              flame_params.rel_vol_[k] = buffer[k];
            }
            if(j==num_states-1){
              y[k*num_states + j] /= flame_params.ref_temperature_;
            }
          }
          // Exit i loop
          break;
        }
        if (i==num_vars_file-1 && state_name != file_state_names[i]) {
          // Variable not found
          printf("WARNING: %s not found in restart file.\n", state_name.c_str() );
        }
      } // for i<num_vars_file
    } // for j<num_states


    // Read mass flux
#ifdef ZERORK_MPI
    disp = 2*sizeof(int) + sizeof(double) + num_vars_file*sizeof(char)*64
      + num_vars_file*sizeof(double) // left BC from prev vars
      + sizeof(double) // left BC from cur var
      + (my_pe + num_vars_file*npes)*num_local_points*sizeof(double);
    MPI_File_set_view(restart_file, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    MPI_File_read(restart_file, &buffer[0], num_local_points, MPI_DOUBLE, MPI_STATUS_IGNORE);
#else
    restart_file.read((char*)&bc, sizeof(double)); //read/skip BC
    restart_file.read((char*)&buffer[0], num_local_points*sizeof(double));
#endif
    for (int k=0; k<num_local_points; ++k) {
      y[k*num_states + num_states-2] = buffer[k];
    }

    //close file
#ifdef ZERORK_MPI
    MPI_File_close(&restart_file);
#else
    restart_file.close();
#endif

  } // if restart_file

  // Set jfix and Tfix
  // TODO: make the parallel implementation cleaner?
  int local_jfix = 1000000;
  for(int j=0; j<num_local_points; ++j) {
    int jglobal = j + my_pe*num_local_points;
    int temp_id = j*num_states+num_species + 1;
    if( y[temp_id] > (flame_params.inlet_temperature_ +
                      flame_params.deltaTfix_/flame_params.ref_temperature_) ) {
      local_jfix = jglobal;
      break;
    }
  }
#if ZERORK_MPI
  MPI_Allreduce(&local_jfix,&flame_params.j_fix_,1,MPI_INT,MPI_MIN,flame_params.comm_);
#else
  flame_params.j_fix_ = local_jfix;
#endif

  double local_Tfix = 0.0;
  for(int j=0; j<num_local_points; ++j) {
    int jglobal = j + my_pe*num_local_points;
    int temp_id = j*num_states+num_species + 1;
    if( jglobal == flame_params.j_fix_ ) local_Tfix = y[temp_id];
  }
#if ZERORK_MPI
  MPI_Allreduce(&local_Tfix,&flame_params.temperature_fix_,1,MPI_DOUBLE,MPI_MAX,flame_params.comm_);
#else
  flame_params.temperature_fix_ = local_Tfix;
#endif

}
