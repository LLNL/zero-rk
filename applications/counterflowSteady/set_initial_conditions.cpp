#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
#include <utilities/string_utilities.h>
#include <utilities/math_utilities.h>
#include <utilities/file_utilities.h>

#include "set_initial_conditions.h"

int sign(const double x)
{
    return (x > 0) ? 1 :
           (x < 0) ? -1
                   : 0;
}

int sign(const int x)
{
    return (x > 0) ? 1 :
           (x < 0) ? -1
                   : 0;
}



// Set initial composition
//void SetInitialComposition(FlameParams &flame_params, double *y, double *ydot, double *id, double *time)
void SetInitialComposition(FlameParams &flame_params, N_Vector yvec, double *time)
{
  const int num_points = flame_params.z_.size();
  const int num_local_points  = flame_params.num_local_points_;
  const int num_states  = flame_params.reactor_->GetNumStates();
  const int num_species = flame_params.reactor_->GetNumSpecies();

  const double pressure = flame_params.parser_->pressure();
  const double ref_temperature = flame_params.parser_->ref_temperature();
  const double ref_momentum = flame_params.parser_->ref_momentum();

  double *y = NV_DATA_P(yvec);

  int my_pe = flame_params.my_pe_;
  MPI_Comm comm = flame_params.comm_;

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

    MPI_File_open(flame_params.comm_, filename2, MPI_MODE_RDONLY, MPI_INFO_NULL, &restart_file);

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
    for (int k=0; k<num_local_points*num_states; ++k)
      y[k] = 0.0;

    for(int j=0; j<num_states; ++j) {
      string state_name = zerork::utilities::GetLowerCase(flame_params.reactor_->GetNameOfStateId(j));
      for(int i=0; i<num_vars_file; ++i) {
        string file_state_name = zerork::utilities::GetLowerCase(file_state_names[i]);
        if(state_name == file_state_name) {
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
            if(j==num_species) {
              flame_params.rel_vol_[k] = buffer[k];
            }
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

    // Read mass flux BC
    disp = 2*sizeof(int) + sizeof(double) + num_vars_file*sizeof(char)*64
      + num_vars_file*2*sizeof(double) // Left & right BCs from previous variables
      + num_vars_file*npes*num_local_points*sizeof(double);
    MPI_File_set_view(restart_file, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    MPI_File_read(restart_file, &buffer[0], 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
    if(!flame_params.parser_->finite_separation())
      flame_params.mass_flux_fuel_ = buffer[0];

    // Read mass flux
    disp = 2*sizeof(int) + sizeof(double) + num_vars_file*sizeof(char)*64
      + num_vars_file*2*sizeof(double) // Left & right BCs from previous variables
      + sizeof(double) // Left BC from current variable
      + (my_pe + num_vars_file*npes)*num_local_points*sizeof(double);
    MPI_File_set_view(restart_file, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    MPI_File_read(restart_file, &buffer[0], num_local_points, MPI_DOUBLE, MPI_STATUS_IGNORE);
    for (int k=0; k<num_local_points; ++k) {
      y[k*num_states + num_species] = buffer[k];
    }

    // Read mass flux BC
    disp = 2*sizeof(int) + sizeof(double) + num_vars_file*sizeof(char)*64
      + num_vars_file*2*sizeof(double) // Left & right BCs from previous variables
      + sizeof(double) // Left BC from current variable
      + (num_vars_file+1)*npes*num_local_points*sizeof(double);
    MPI_File_set_view(restart_file, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
    MPI_File_read(restart_file, &buffer[0], 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
    if(!flame_params.parser_->finite_separation())
      flame_params.mass_flux_oxidizer_ = buffer[0];

    MPI_File_close(&restart_file);

  } else {
    printf("# Steady solver requires restart_file to be specified.\n");
    exit(-1);
  }// if restart_file

  // Recompute relative volume?
  const double RuTref_p = flame_params.reactor_->GetGasConstant()*
    flame_params.ref_temperature_/flame_params.parser_->pressure();

  for(int j=0; j<num_local_points; ++j) {
    int temp_id = j*num_states+num_species + 1;

    double mass_fraction_sum = 0.0;
    for(int k=0; k<num_species; ++k)
      mass_fraction_sum += flame_params.inv_molecular_mass_[k]*y[j*num_states+k];

    flame_params.rel_vol_[j] = RuTref_p*y[temp_id]*mass_fraction_sum;
  }

  // Compute stagnation plane location
  // ONLY WORKS IN SERIAL FOR NOW
  long int dsize;
  int jContBC;
  double *mbuf, *mloc;
  mloc = (double *)malloc(num_local_points*sizeof(double));
  if(my_pe == 0)
    mbuf = (double *)malloc(num_points*sizeof(double));

  for(int j=0; j<num_local_points; j++)
    mloc[j]= y[j*num_states + num_species];

  // Gather global mass flux on root
  dsize = num_local_points;
  MPI_Gather(mloc,
             dsize,
             PVEC_REAL_MPI_TYPE,
             mbuf,
             dsize,
             PVEC_REAL_MPI_TYPE,
             0,
             comm);

  if(my_pe == 0) {
    if(flame_params.flame_type_ == 0 || flame_params.flame_type_ == 2) {
      // Find grid point closest to old stagnation point location
      int jStart;
      for(int j=0; j<num_points; j++) {
        if(flame_params.z_[j] > flame_params.stagnation_plane_) {
          jStart = j;
          break;
        }
      }
      jContBC = jStart;
      for(int i=1; jStart+i < num_points || jStart >= i; i++) {
        if(jStart+i < num_points && sign(mbuf[jStart+i]) != sign(mbuf[jStart])) {
          jContBC = jStart + i;
          break;
        } else if (jStart>=i && sign(mbuf[jStart-i]) != sign(mbuf[jStart])) {
          jContBC = jStart-i+1;
          break;
        }
      }
      if(jContBC == 0) {
        assert(mbuf[jContBC] <=0);
        flame_params.stagnation_plane_ = flame_params.z_[0] - mbuf[0]*
          flame_params.dz_[0]/(mbuf[1]-mbuf[0]);
      } else if (jContBC == num_points-1) {
        assert(mbuf[jContBC] >=0);
        flame_params.stagnation_plane_ = flame_params.z_[jContBC] - mbuf[jContBC]*
        flame_params.dz_[jContBC]/(mbuf[jContBC]-mbuf[jContBC-1]);
      } else {
        assert(mbuf[jContBC]*mbuf[jContBC-1] <=0); // test opposite sign
        flame_params.stagnation_plane_ = flame_params.z_[jContBC] - mbuf[jContBC]*
          flame_params.dz_[jContBC]/(mbuf[jContBC]-mbuf[jContBC-1]);
      }
      jContBC -= 1; //??
    } else if (flame_params.flame_type_ == 1) {
      flame_params.stagnation_plane_ = flame_params.length_;
      jContBC = flame_params.num_points_-1;
    }
  }
  // Broadcast stagnation plane location and index
  MPI_Bcast(&flame_params.stagnation_plane_, 1, MPI_DOUBLE, 0, comm);
  MPI_Bcast(&jContBC, 1, MPI_INT, 0, comm);
  flame_params.jContBC_ = jContBC;
  if(my_pe == 0)
    printf("# Stagnation plane at x: %5.3e, j: %d\n", flame_params.stagnation_plane_,
           flame_params.jContBC_);

  // Set fixed temperature to initial temperature
  for(int j=0; j<num_local_points; j++) {
    flame_params.fixed_temperature_[j] = y[j*num_states + num_species + 1];
  }

}
