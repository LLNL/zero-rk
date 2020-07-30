#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <string>
#include <vector>
#include <map>

#include "zerork/elemental_composition.h"

#include "janaf_thermo.h"
#include "thermo_parser.h"
#include "thermo_fix.h"
#include "ThermoCheckIFP.h"

int main(int argc, char *argv[])
{
  // ------------------------------------------
  std::vector<std::string> species_name;
  std::vector<std::string> mech_species_name;
  std::vector<std::string> species_definition;
  std::vector<zerork::ElementalComposition> species_composition;
  std::vector<JanafThermoData> janaf_data;
  std::vector<double> global_T_range; 

  std::map<std::string, int> species_position;
  std::vector<int> repair_species_id;

  int num_extrema, num_extrema_species;
  std::vector<double> Cp_extrema;
  std::vector<double> T_extrema;
  int species_id_max_jump[3];
  int num_thermo_species, num_mech_species;
  
  if(argc != 2) {
    printf("ERROR: incorrect command line usage.\n");
    printf("       use instead %s <input file parameters>\n",argv[0]);
    exit(-1);
  }
  ThermoCheckIFP parser(argv[1]);

  num_thermo_species = ParseThermoFile(parser.therm_file().c_str(),
                                       species_name,
                                       species_definition,
                                       species_composition,
                                       janaf_data,
                                       global_T_range);
  if(num_thermo_species == 0) {
    printf("ERROR: ParseThermoFile(...) found zero species\n");
    exit(-1);
  }
  printf("! Found %d distinct species in thermo file %s\n",
         num_thermo_species,parser.therm_file().c_str());

  num_mech_species = 0;
  if(parser.use_min_species()) {
    num_mech_species = ParseMechFileForSpecies(parser.mech_file().c_str(),
                                               mech_species_name);
    printf("! Found %d distinct species in mechanism file %s\n",
           num_mech_species,parser.mech_file().c_str());
  }
  
  if(num_mech_species > 0 && num_mech_species <= num_thermo_species) {
    // construct the species_position map for the species found in the thermo
    // file
    species_position.clear();
 
    for(int j=0; j<num_mech_species; ++j) {
      species_position[mech_species_name[j]] = j;
    }
    for(int j=0; j<num_thermo_species; ++j) {
      if(species_position.find(species_name[j]) == species_position.end()) {
        // thermodynamic species j is not found in the mechanism file
        printf("! INFO: ignoring species %s from the thermo file %s\n",
               species_name[j].c_str(),
               parser.therm_file().c_str());
        printf("        that is not found in the mechanism file %s\n",
               parser.mech_file().c_str()); 
      } else {
        // thermodynamic species j is found in the mechanism file
        repair_species_id.push_back(j);
      }
    }
  } else if(num_mech_species > num_thermo_species) {

    printf("INFO: number of species %d in the mechanism file %s\n",
           num_mech_species,parser.mech_file().c_str());
    printf("      exceeds the number of species %d in the thermo file %s\n",
           num_thermo_species,parser.therm_file().c_str());
    printf("      Repairing all %d thermodynamic species\n",
           num_thermo_species);
    fflush(stdout);
    for(int j=0; j<num_thermo_species; ++j) {
      repair_species_id.push_back(j);
    } 
  }
  else {
    // no mechanism species found
    printf("INFO: no species found in mechanism file %s\n",
           parser.mech_file().c_str());
    printf("      Repairing all %d thermodynamic species\n",
           num_thermo_species);
    fflush(stdout);
    for(int j=0; j<num_thermo_species; ++j) {
      repair_species_id.push_back(j);
    }
  }

  printf("Scanning thermodynamics for local extrema in Cp/R\n");
  fflush(stdout);
  num_extrema_species = 0;
  for(size_t j=0; j<repair_species_id.size(); ++j) {

    int spec_id = repair_species_id[j];
    // The Janaf thermodynamics will be reported as having an extrema
    // if:
    //
    num_extrema = GetJanafExtrema(janaf_data[spec_id],
                                  1.0e-12, // imaginary root relative tol
                                  1.0e-16, // imaginary root absolute tol
                                  &T_extrema,
                                  &Cp_extrema);
    if(num_extrema > 0) {
      printf("! Extrema found for species %d [%s]\n",
             spec_id,
             species_name[spec_id].c_str());
      for(int k=0; k<num_extrema; ++k) {
        printf("!   at T = %10.5f, Cp/R = %10.5f\n",
               T_extrema[k],Cp_extrema[k]);
        fflush(stdout);
      }
      ++num_extrema_species;
    }
  }
  printf("A total of %d species thermodynamics are found\n",
         num_extrema_species);
  printf("to be non-monotonic\n");

  return 0;
}
