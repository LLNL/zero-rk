#include <stdio.h>

#include <vector>
#include <string>
#include <cassert>

#include <utilities/file_utilities.h>
#include <utilities/string_utilities.h>

#ifdef ZERORK_MPI
#include <mpi.h>
#endif

#include "zerork_flame_api.h"
#include "ZeroRKFlameAPITesterIFP.h"

#include "zerork/mechanism.h"

static zerork_flame_handle zf_handle;

void check_status(zerork_flame_status_t status) {
  if(status != ZERORK_FLAME_STATUS_SUCCESS) {
    printf("flame_api returned error : %d\n", status);
    exit(2);
  }
}

int ReadProfileFile(const char* file, int num_species, std::vector<double>& grid, std::vector<double>& T, std::vector<double>& mass_fractions);
int WriteProfileFile(const char* file, int num_species, std::vector<double>& grid, std::vector<double>& T, std::vector<double>& mass_fractions);

int main(int argc, char *argv[])
{
  int nranks = 1;
  int myrank = 0;
#ifdef ZERORK_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

  if(argc < 2) {
    printf("# ERROR: Incorrect command line usage.\n");
    printf("#        use %s <input parameters>\n",argv[0]);
    exit(1);
  }

  ZeroRKFlameAPITesterIFP inputFileDB(argv[1]);

  const char* mechfilename = inputFileDB.mechanism_file().c_str();
  const char* thermfilename = inputFileDB.therm_file().c_str();
  const char* transfilename = inputFileDB.transport_file().c_str();

  zf_handle = zerork_flame_init();

  check_status(zerork_flame_read_options_file(argv[1], zf_handle));
  check_status(zerork_flame_set_int_option("verbosity", inputFileDB.verbosity(), zf_handle));
  check_status(zerork_flame_set_double_option("relative_tolerance", inputFileDB.relative_tolerance(), zf_handle));
  check_status(zerork_flame_set_double_option("absolute_tolerance", inputFileDB.absolute_tolerance(), zf_handle));
  check_status(zerork_flame_set_int_option("pseudo_unsteady", inputFileDB.pseudo_unsteady(), zf_handle));
  check_status(zerork_flame_set_double_option("pseudo_unsteady_dt", inputFileDB.pseudo_unsteady_dt(), zf_handle));
  check_status(zerork_flame_set_double_option("pseudo_unsteady_min_dt", inputFileDB.pseudo_unsteady_min_dt(), zf_handle));
  check_status(zerork_flame_set_double_option("pseudo_unsteady_max_dt", inputFileDB.pseudo_unsteady_max_dt(), zf_handle));
  check_status(zerork_flame_set_int_option("integrator_type", inputFileDB.integrator_type(), zf_handle));
  check_status(zerork_flame_set_int_option("convective_scheme_type", inputFileDB.convective_scheme_type(), zf_handle));
  check_status(zerork_flame_set_string_option("transport_model", inputFileDB.transport_model().c_str(), zf_handle));
  check_status(zerork_flame_set_string_option("constant_lewis_setting", inputFileDB.constant_lewis_setting().c_str(), zf_handle));
  check_status(zerork_flame_set_int_option("constant_lewis_grid_point", inputFileDB.constant_lewis_grid_point(), zf_handle));
  check_status(zerork_flame_set_string_option("temperature_fix_setting", inputFileDB.temperature_fix_setting().c_str(), zf_handle));
  check_status(zerork_flame_set_double_option("temperature_fix_value", inputFileDB.temperature_fix_value(), zf_handle));
  check_status(zerork_flame_set_string_option("kinsol_strategy", inputFileDB.kinsol_strategy().c_str(), zf_handle));
  check_status(zerork_flame_set_input_files(mechfilename, thermfilename, transfilename, zf_handle));
  check_status(zerork_flame_load_mechanism(zf_handle));

  const char* cklogfilename = zerork::utilities::null_filename; //We already parsed in load_mechanism, no need to have another log
  zerork::mechanism mech(mechfilename, thermfilename, cklogfilename);
  int num_species = mech.getNumSpecies();

  double flame_speed = inputFileDB.flame_speed();
  double pressure = inputFileDB.pressure();
  const char* flame_start_profile = inputFileDB.flame_start_profile().c_str();
  std::vector<double> grid;
  std::vector<double> T;
  std::vector<double> mass_fractions;
  int num_grid_points = 0;
  if(myrank == 0) {
    int read_status = ReadProfileFile(flame_start_profile, num_species, grid, T, mass_fractions);
    if(read_status != 0) {
      printf("Error reading flame_start_profile (%s). Exiting\n",flame_start_profile);
      return 1;
    }
    num_grid_points = grid.size();
    assert(T.size() == num_grid_points);
    assert(mass_fractions.size() == num_grid_points*num_species);
  }
  
  //GT: We can separate out the grid-info and pressure from this call, but decided this was simpler
  check_status(zerork_flame_solve(num_grid_points,       //scalar (in)
                                  grid.data(),           //array (in)
                                  pressure,              //scalar (in)
                                  &flame_speed,          //scalar (in/out)
                                  T.data(),              //array (in/out) of temperatures at each grid point
                                  mass_fractions.data(), //array (in/out) of mass fractions ([grid_idx*num_species + species_idx])
                                  zf_handle));

  check_status(zerork_flame_free(zf_handle));

  if(myrank == 0) {
    printf("final flame speed: %g\n",flame_speed);

    const char* flame_end_profile = inputFileDB.flame_end_profile().c_str();
    WriteProfileFile(flame_end_profile, num_species, grid, T, mass_fractions);
  }

#ifdef ZERORK_MPI
  MPI_Finalize();
#endif
}

int ReadProfileFile(const char* file, int num_species, std::vector<double>& grid, std::vector<double>& T, std::vector<double>& mass_fractions) {
  std::ifstream profile_file;
  std::string line;
  std::string delimiters = std::string(zerork::utilities::WHITESPACE) + ",";
  std::vector<std::string> fields;

  profile_file.open(file);

  if(profile_file) {
    while(zerork::utilities::GetAnyLine(profile_file,&line)) {
      zerork::utilities::SplitStringToVector(line, delimiters, &fields);
      if(fields.size() == 2+num_species) {
        grid.push_back(atof(fields[0].c_str()));
        for(int k=0; k<num_species; ++k) {
          mass_fractions.push_back(atof(fields[1+k].c_str()));
        }
        T.push_back(atof(fields[num_species+1].c_str()));
      }
    }
  } else {
    return 1;
  }
  return 0;
}

int WriteProfileFile(const char* file, int num_species, std::vector<double>& grid, std::vector<double>& T, std::vector<double>& mass_fractions) {
  const int num_grid_points = grid.size();
  assert(num_grid_points == T.size());
  assert(num_grid_points*num_species == mass_fractions.size());

  std::ofstream profile_file;
  profile_file.open(file);

  if(profile_file) {
    for(int j = 0; j < num_grid_points; ++j) {
      profile_file << grid[j];
      for(int k = 0; k < num_species; ++k) {
        profile_file << "\t" << mass_fractions[j*num_species+k];
      }
      profile_file << "\t" << T[j];
      profile_file << std::endl;
    }
  }
  return 0;
}

