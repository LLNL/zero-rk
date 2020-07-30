#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "zerork/elemental_composition.h"

#include "janaf_thermo.h"
#include "thermo_parser.h"
#include "thermo_fix.h"
#include "ThermoCheckIFP.h"

void WriteJanafSpeciesToFile(FILE *fptr,
                             std::string &species_definition,
                             JanafThermoData &data);
void WriteThermoComparisonToFile(FILE *fptr,
                                 int record_num,
                                 int num_print_points,
                                 std::string &species_name,
                                 double Cp_R_max_diff,
                                 double H_RT_max_diff,
                                 double S_R_max_diff,
                                 JanafThermoData &data1,
                                 JanafThermoData &data2);

std::vector<std::string> GetElemVector(std::vector<zerork::ElementalComposition>& species_composition);

void WriteThermoSummaryHeader(FILE *fptr,
                             const std::vector<std::string> &elems);
void WriteThermoSummarySpecies(FILE *fptr,
                              const std::vector<std::string> &elems,
                              const std::string &species_name,
                              const zerork::ElementalComposition &species_composition,
                              const JanafThermoData &data);

void WriteIsomerJanafBlock(FILE* fptr,
                           int num_print_points,
                           const zerork::ElementalComposition &composition,
                           const std::vector<std::string> &isomer_species_names,
                           const std::vector<JanafThermoData> &isomer_janaf_data);

double MaxDouble(const double a, const double b) {
  return ((a>b) ? a : b);
}

typedef struct {
  double specific_heat;
  double enthalpy;
  double entropy;
  int species_id;
} MaxThermoChange;

double AbsMaxThermoChange(const MaxThermoChange *aptr) {
  return MaxDouble(fabs(aptr->specific_heat),
                        MaxDouble(fabs(aptr->enthalpy),
                                  fabs(aptr->entropy)));
}

int CompareMaxThermoChange(const void *a, const void *b) {
  MaxThermoChange *aptr = (MaxThermoChange *) a;
  MaxThermoChange *bptr = (MaxThermoChange *) b;
 
  double amax = AbsMaxThermoChange(aptr);
  double bmax = AbsMaxThermoChange(bptr);
  if(amax > bmax) {
    return -1; // 'a' goes first if larger
  } else if(amax < bmax) {
    return 1;  // 'b' goes first if larger
  }
  return 0;
}


int main(int argc, char *argv[])
{
  // TODO: parameters to be specified in input file

  // ------------------------------------------
  double T_fixed, thermo_prob_atol;
  int num_prob_thermo=0;
  FILE *refit_fptr,*thermo_prob_fptr;  
  std::vector<std::string> species_name;
  std::vector<std::string> mech_species_name;
  std::vector<std::string> species_definition;
  std::vector<zerork::ElementalComposition> species_composition;
  std::vector<JanafThermoData> janaf_data;
  std::vector<double> global_T_range; 

  std::map<int, int> repair_species_map;
  std::map<std::string, int> species_position;
  std::map<double, int> max_thermo_difference;
  std::vector<int> repair_species_id;
  MaxThermoChange *repair_diff;

  std::vector<JanafThermoData> refit_janaf;
  int num_extrema, num_extrema_species;
  std::vector<double> Cp_extrema;
  std::vector<double> T_extrema;

  double Cp_R_max_diff,TatCp_R_max_diff;
  double H_RT_max_diff,TatH_RT_max_diff;
  double S_R_max_diff, TatS_R_max_diff;
  double T_match_jump[3];
  double max_jump[3];
  int species_id_max_jump[3];
  int num_thermo_species, num_mech_species;
  bool is_species_uppercase = false;
  bool try_refit_T_match = false;
  
  if(argc != 2) {
    printf("ERROR: incorrect command line usage.\n");
    printf("       use instead %s <input file parameters>\n",argv[0]);
    exit(-1);
  }
  ThermoCheckIFP parser(argv[1]);

  refit_fptr = fopen(parser.repair_file().c_str(),"w");
  if(refit_fptr == NULL) {
    printf("ERROR: could not open file %s for write operation\n",
           parser.repair_file().c_str());
    exit(-1);
  }
  thermo_prob_fptr = fopen(parser.change_file().c_str(),"w");
  if(thermo_prob_fptr == NULL) {
    printf("ERROR: could not open file %s for write operation\n",
           parser.change_file().c_str());
    fclose(refit_fptr);
    exit(-1);
  }
  thermo_prob_atol = parser.repair_atol();
  T_fixed = parser.retain_temperature();
  try_refit_T_match = parser.find_best_Tmatch();
  
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
  if(num_thermo_species == -1) {
    printf("ERROR: ParseThermoFile(...) failed\n");
    exit(-1);
  }
  printf("! Found %d distinct species in thermo file %s\n",
         num_thermo_species,parser.therm_file().c_str());
  // copy original janaf data to refit_data
  refit_janaf.clear();
  for(int j=0; j<num_thermo_species; ++j) {
    refit_janaf.push_back(janaf_data[j]);
    PrintJanaf(janaf_data[j]);
  }

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
        repair_species_map[j] = repair_species_id.back();
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
      repair_species_map[j] = repair_species_id.back();
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
      repair_species_map[j] = repair_species_id.back();
    }
  }

  // this won't work for species like '2ma'
  if('A' <= species_name[0][0] && species_name[0][0] <= 'Z') {
    is_species_uppercase = true;
  }

  printf("Global T range [%6.2f, %6.2f, %6.2f]\n",
         global_T_range[0],
         global_T_range[1],
         global_T_range[2]);
  if(is_species_uppercase) {
    fprintf(refit_fptr,"THERMO\n");
  } else {
    fprintf(refit_fptr,"thermo\n");
  }
  // Global Temperature range has fixed width of 10 columns
  // Note that the match temperature is restricted to 2 decimal places, 
  // which is retained here
  fprintf(refit_fptr,"%10.4f%10.2f%10.4f\n",
          global_T_range[0],
          global_T_range[1],
          global_T_range[2]);

  // initialize the jump and thermo difference maxima
  for(int j=0; j<3; ++j) {
    max_jump[j]=-1.0e300;
  }

  repair_diff = new MaxThermoChange[repair_species_id.size()];

  for(size_t k=0; k<repair_species_id.size(); ++k) {
    int j = repair_species_id[k];
    printf("! Refitting species %s:\n",species_name[j].c_str());
    if(parser.add_comments()) {
      fprintf(refit_fptr,
              "! Refitting species %s:\n",
              species_name[j].c_str());
    }
    if(parser.fix_Tmatch()) {
      // refit the lower temperature polynomial with the globally defined
      // matching temperature specified on line 2 of the thermo file
      RefitKeepHighGlobalTMatch(parser.num_repair_points(),
                                2, // number of decimal places in new temp
                                T_fixed,
                                &janaf_data[j],
                                global_T_range[1], // use global T_match
                                &refit_janaf[j]);
    } else {
      // refit the lower temperature polynomial using one of two methods:
      // try_refit_T_match = false : use the same matching temperature as the 
      //                             original species definition
      // try_refit_T_match = true  : search for a matching temperature that 
      //                             minimizes the Cp/R jump at T_match
      RefitKeepHigh(parser.num_repair_points(),
                  2,                // number of decimal places in new temp
                  T_fixed,
                  &janaf_data[j],
                  try_refit_T_match,  // recompute better matching temperature
                  &refit_janaf[j]);
    }
    // compute Cp difference
    Cp_R_max_diff = MaxDifference(JanafSpecificHeat,
                                  parser.num_test_points(),
                                  &janaf_data[j],
                                  &refit_janaf[j],
                                  &TatCp_R_max_diff);

    printf("!   max difference Cp/R = %10.4e at T = %10.4f K\n",
           Cp_R_max_diff,TatCp_R_max_diff);
    if(parser.add_comments()) {
      fprintf(refit_fptr,
              "!   max difference Cp/R = %10.4e at T = %10.4f K\n",
              Cp_R_max_diff,TatCp_R_max_diff);
    }
    // compute enthalpy difference
    H_RT_max_diff = MaxDifference(JanafEnthalpy,
                                  parser.num_test_points(),
                                  &janaf_data[j],
                                  &refit_janaf[j],
                                  &TatH_RT_max_diff);
    printf("!   max difference H/RT = %10.4e at T = %10.4f K\n",
           H_RT_max_diff,TatH_RT_max_diff);
    if(parser.add_comments()) {
      fprintf(refit_fptr,
              "!   max difference H/RT = %10.4e at T = %10.4f K\n",
              H_RT_max_diff,TatH_RT_max_diff);
    }
    // compute entropy difference
    S_R_max_diff = MaxDifference(JanafEntropy,
                                 parser.num_test_points(),
                                 &janaf_data[j],
                                 &refit_janaf[j],
                                 &TatS_R_max_diff);

    printf("!   max difference S/R  = %10.4e at T = %10.4f K\n",
           S_R_max_diff,TatS_R_max_diff);
    if(parser.add_comments()) {
      fprintf(refit_fptr,
              "!   max difference S/R  = %10.4e at T = %10.4f K\n",
              S_R_max_diff,TatS_R_max_diff);
    }

    repair_diff[k].specific_heat = Cp_R_max_diff;
    repair_diff[k].enthalpy      = H_RT_max_diff;
    repair_diff[k].entropy       = S_R_max_diff;
    repair_diff[k].species_id    = j;

    // compute the jump (high - low)
    JanafJump(refit_janaf[j],T_match_jump);
    printf("!   T_match = %12.6f jump\n",refit_janaf[j].T_match);
    printf("!   Delta(Cp/R, H/RT, S/R) = (%10.3e,%10.3e,%10.3e)\n",
           T_match_jump[0],T_match_jump[1],T_match_jump[2]);
 
    if(parser.add_comments()) {
      fprintf(refit_fptr,
              "!   T_match = %12.6f jump\n",refit_janaf[j].T_match);
      fprintf(refit_fptr,
              "!   Delta(Cp/R, H/RT, S/R) = (%10.3e,%10.3e,%10.3e)\n",
              T_match_jump[0],T_match_jump[1],T_match_jump[2]);
    }

    for(int m=0; m<3; ++m) {
      if(fabs(T_match_jump[m]) > max_jump[m]) {
        max_jump[m] = fabs(T_match_jump[m]);
        species_id_max_jump[m] = j;
      }
    }

    if(AbsMaxThermoChange(&repair_diff[k]) > thermo_prob_atol) {
      ++num_prob_thermo;
    }

    WriteJanafSpeciesToFile(refit_fptr,
                            species_definition[j],
                            refit_janaf[j]);
    // printf("%s: ", species_name[j].c_str()); PrintJanaf(janaf_data[j]); fflush(stdout);
  }
  if(is_species_uppercase) {
    fprintf(refit_fptr,"END\n");
  } else {
    fprintf(refit_fptr,"end\n");
  }
  // finished writing the repaired thermodynamics file
  fclose(refit_fptr);

  qsort(&repair_diff[0],
        repair_species_id.size(),
        sizeof(MaxThermoChange),
        CompareMaxThermoChange);

  if(num_prob_thermo > 0) {
    // record the information for the problem thermo file
    fprintf(thermo_prob_fptr,
            "# JANAF 1 - original thermo definition read from: %s\n",
            parser.therm_file().c_str());
    fprintf(thermo_prob_fptr,
            "# JANAF 2 - refit thermo definitio written to   : %s\n",
            parser.repair_file().c_str());
    fprintf(thermo_prob_fptr,
            "#           refit T_fixed   [K] : %.18g\n",T_fixed);
    fprintf(thermo_prob_fptr,
            "#           number of fit points: %d\n",
            parser.num_repair_points());
  }
  for(int j=0; j<num_prob_thermo; ++j) {
    
    int species_id = repair_diff[j].species_id;

    WriteThermoComparisonToFile(thermo_prob_fptr,
                                j, //record_num
                                parser.num_print_points(),
                                species_name[species_id],
                                repair_diff[j].specific_heat,
                                repair_diff[j].enthalpy,
                                repair_diff[j].entropy,
                                janaf_data[species_id],
                                refit_janaf[species_id]);
  }
  // report the extrema
  printf("\n\nThermo Refit Extrema:\n");
  printf("Largest absolute change in the repaired thermodynamics\n");
  printf("  species %d [%s]\n",
         repair_diff[0].species_id,
         species_name[repair_diff[0].species_id].c_str());
  printf("  Cp/R difference   = %10.4e\n",
         repair_diff[0].specific_heat);
  printf("  H/RT difference   = %10.4e\n",
         repair_diff[0].enthalpy);

  printf("  S/R  difference   = %10.4e\n",
         repair_diff[0].entropy);

  printf("  max Cp/R T_match jump = %10.4e for species %d [%s]\n",
         max_jump[0],
         species_id_max_jump[0],
         species_name[species_id_max_jump[0]].c_str());
  printf("  max H/RT T_match jump = %10.4e for species %d [%s]\n",
         max_jump[1],
         species_id_max_jump[1],
         species_name[species_id_max_jump[1]].c_str());
  printf("  max S/R  T_match jump = %10.4e for species %d [%s]\n",
         max_jump[2],
         species_id_max_jump[2],
         species_name[species_id_max_jump[2]].c_str());

  printf("%d thermo fits with a difference greater than atol = %8.2e\n",
         num_prob_thermo,parser.repair_atol());
  printf("written to file: %s\n",parser.change_file().c_str());

  printf("Scanning refit thermodynamics for local extrema in Cp/R\n");
  fflush(stdout);
  num_extrema_species = 0;
  for(size_t j=0; j<repair_species_id.size(); ++j) {
    int spec_id = repair_species_id[j];
    // The Janaf thermodynamics will be reported as having an extrema
    // if:
    //
    num_extrema = GetJanafExtrema(refit_janaf[spec_id],
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

      printf("! Writing Non-monotonic species to problem file.\n");
      WriteThermoComparisonToFile(thermo_prob_fptr,
                                  num_prob_thermo+num_extrema_species, //record_num
                                  parser.num_print_points(),
                                  species_name[spec_id],
                                  repair_diff[j].specific_heat,
                                  repair_diff[j].enthalpy,
                                  repair_diff[j].entropy,
                                  janaf_data[spec_id],
                                  refit_janaf[spec_id]);
    }
  }
  printf("A total of %d refit species thermodynamics are found\n",
         num_extrema_species);
  printf("to be non-monotonic\n");
  fclose(thermo_prob_fptr);

  //Detailed table of reference enthalpy, entropy and Cp's
  if(parser.detailed_table_file() != std::string("")) {
    FILE* detailed_table_fptr = fopen(parser.detailed_table_file().c_str(),"w");
    if(detailed_table_fptr == NULL) {
      printf("ERROR: could not open file %s for write operation\n",
             parser.detailed_table_file().c_str());
      exit(-1);
    }

    typedef std::map<zerork::ElementalComposition,size_t> comp_map_idx_t;
    comp_map_idx_t  species_composition_idx_map;
    for(size_t i=0; i<species_composition.size(); ++i)
    {
      species_composition_idx_map[species_composition[i]]=i;
    }

    std::vector<std::string> elems = GetElemVector(species_composition);
    WriteThermoSummaryHeader(detailed_table_fptr, elems);
    

    for(comp_map_idx_t::const_iterator itr=species_composition_idx_map.begin();
        itr != species_composition_idx_map.end(); ++itr) {
      size_t species_id = itr->second;
      if(repair_species_map.find(species_id) == repair_species_map.end()) {
        continue;
      }

      WriteThermoSummarySpecies(detailed_table_fptr,
                         elems,
                         species_name[species_id],
                         itr->first,
                         refit_janaf[species_id]);
    
    }
    fclose(detailed_table_fptr);
  }

  //Isomer plot data
  if(parser.isomer_plot_file() != std::string("")) {
    FILE* isomer_plot_fptr = fopen(parser.isomer_plot_file().c_str(),"w");
    if(isomer_plot_fptr == NULL) {
      printf("ERROR: could not open file %s for write operation\n",
             parser.isomer_plot_file().c_str());
      exit(-1);
    }

    typedef std::map<zerork::ElementalComposition,size_t> comp_map_idx_t;
    comp_map_idx_t species_composition_idx_map;
    for(size_t i=0; i<species_composition.size(); ++i)
    {
      species_composition_idx_map[species_composition[i]]=i;
    }

    std::vector<std::string> elems = GetElemVector(species_composition);

    std::map<std::string,std::vector<int> > composition_map;
    for(comp_map_idx_t::const_iterator itr=species_composition_idx_map.begin();
        itr != species_composition_idx_map.end(); ++itr) {
      size_t species_id = itr->second;
      if(repair_species_map.find(species_id) == repair_species_map.end()) {
        continue;
      }
      std::string comp_str = itr->first.ToString();
      composition_map[comp_str].push_back(species_id);
    }

    for(comp_map_idx_t::const_iterator itr=species_composition_idx_map.begin();
        itr != species_composition_idx_map.end(); ++itr) {
      size_t species_id = itr->second;
      if(repair_species_map.find(species_id) == repair_species_map.end()) {
        continue;
      }
      std::string comp_str = itr->first.ToString();
      if(composition_map[comp_str].size() > 0) {
        std::vector<std::string> isomer_species_names;
        std::vector<JanafThermoData> isomer_janaf_data;
        for(size_t k = 0; k < composition_map[comp_str].size(); ++k) {
          int species_id = composition_map[comp_str][k];
          isomer_species_names.push_back(species_name[species_id]);
          isomer_janaf_data.push_back(refit_janaf[species_id]);
        }
        WriteIsomerJanafBlock(isomer_plot_fptr,
                              parser.num_print_points(),
                              itr->first,
                              isomer_species_names,
                              isomer_janaf_data);
        composition_map[comp_str].clear();
      }
    }
    
    fclose(isomer_plot_fptr);
  }

  delete [] repair_diff;
  return 0;
}

void WriteJanafSpeciesToFile(FILE *fptr,
                             std::string &species_definition,
                             JanafThermoData &data)
{
  size_t definition_length = species_definition.size();
  if(definition_length != 45) {
    // TODO correct the problem
    printf("ERROR: species definition length %lu != 45\n",
           definition_length);
    exit(-1);
  }
  // species line - record 1
  fprintf(fptr,"%s %9.4f %9.4f %7.2f %6d\n",
          species_definition.c_str(),
          data.T_min,
          data.T_max,
          data.T_match,
          1); // thermo line number
  // thermo data - record 2
  fprintf(fptr,"%15.8e%15.8e%15.8e%15.8e%15.8e %4d\n",
          data.high_coef[0],
          data.high_coef[1],
          data.high_coef[2],
          data.high_coef[3],
          data.high_coef[4],
          2); // thermo line number
  // thermo data - record 3
  fprintf(fptr,"%15.8e%15.8e%15.8e%15.8e%15.8e %4d\n",
          data.high_coef[5],
          data.high_coef[6],
          data.low_coef[0],
          data.low_coef[1],
          data.low_coef[2],
          3); // thermo line number
  // thermo data - record 4
  fprintf(fptr,"%15.8e%15.8e%15.8e%15.8e %19d\n",
          data.low_coef[3],
          data.low_coef[4],
          data.low_coef[5],
          data.low_coef[6],
          4); // thermo line number
}
void WriteThermoComparisonToFile(FILE *fptr,
                                 int record_num,
                                 int num_print_points,
                                 std::string &species_name,
                                 double Cp_R_max_diff,
                                 double H_RT_max_diff,
                                 double S_R_max_diff,
                                 JanafThermoData &data1,
                                 JanafThermoData &data2)
{
  const double T_epsilon = 1.0e-8;
  double T_min = data1.T_min;
  double T_max = data1.T_max;
  double T_match = data1.T_match;
  double T_plot;

  if(data2.T_min < T_min) {
    T_min = data2.T_min;
  }
  if(data2.T_max > T_max) {
    T_max = data2.T_max;
  }

  fprintf(fptr,"#-----------------------------------------------------------------------------\n");
  fprintf(fptr,"# Block Number: %d\n",record_num);
  fprintf(fptr,
          "# Species: %s T_match (1): %7.2f (2): %7.2f\n#\n",
          species_name.c_str(),
          data1.T_match,
          data2.T_match);
  fprintf(fptr,
          "#   max Cp/R difference = %12.4e\n",
          Cp_R_max_diff);
  fprintf(fptr,
          "#   max H/RT difference = %12.4e\n",
          H_RT_max_diff);
  fprintf(fptr,
          "#   max S/R  difference = %12.4e\n#\n",
          S_R_max_diff);
  fprintf(fptr,
          "#   Temp [K]      Cp/R JANAF 1      H/RT JANAF 1      S/R  JANAF 1      Cp/R JANAF 2      H/RT JANAF 2      S/R  JANAF 2\n");

  for(int j=0; j<num_print_points; ++j) {
    T_plot = T_min + 
      (T_match - T_min)*static_cast<double>(j)/(num_print_points-1.0)
      - T_epsilon; // if(T_plot == T_match), then higher branch is used
                   // T_epsilon ensures that a point close to the matching
                   // temperature is used for the lower branch
    fprintf(fptr,
            "%12.5f  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e\n",
            T_plot,
            JanafSpecificHeat(T_plot,&data1),
            JanafEnthalpy(T_plot,&data1),
            JanafEntropy(T_plot,&data1),
            JanafSpecificHeat(T_plot,&data2),
            JanafEnthalpy(T_plot,&data2),
            JanafEntropy(T_plot,&data2));    
  }
  for(int j=0; j<num_print_points; ++j) {
    T_plot = T_match + 
      (T_max - T_match)*static_cast<double>(j)/(num_print_points-1.0);
    fprintf(fptr,
            "%12.5f  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e  %16.7e\n",
            T_plot,
            JanafSpecificHeat(T_plot,&data1),
            JanafEnthalpy(T_plot,&data1),
            JanafEntropy(T_plot,&data1),
            JanafSpecificHeat(T_plot,&data2),
            JanafEnthalpy(T_plot,&data2),
            JanafEntropy(T_plot,&data2));    
  }
  fprintf(fptr,"\n\n");
}

std::vector<std::string> GetElemVector(std::vector<zerork::ElementalComposition>& species_composition)
{
  zerork::ElementalComposition global_comp;
  for(size_t i = 0; i < species_composition.size(); ++i) {
    global_comp += species_composition[i];
  }
  std::vector<std::string> elems = global_comp.GetElementVector();
  return elems;
}

void WriteThermoSummaryHeader(FILE *fptr,
                             const std::vector<std::string> &elems)
{
  fprintf(fptr,"%s\t","Species");
  for(size_t i = 0; i < elems.size(); ++i)
  {
    fprintf(fptr,"%s\t",elems[i].c_str());
  }
  fprintf(fptr,"H(298)\tS(298)\tCp(300)\tCp(400)\tCp(500)\tCp(600)\tCp(700)\tCp(800)\tCp(1000)\tCp(1500)\n");
}

void WriteThermoSummarySpecies(FILE *fptr,
                              const std::vector<std::string> &elems,
                              const std::string &species_name,
                              const zerork::ElementalComposition &species_composition,
                              const JanafThermoData &data)
{
  fprintf(fptr,"%s\t",species_name.c_str());
  for(size_t i = 0; i < elems.size(); ++i)
  {
    fprintf(fptr,"%d\t",species_composition[elems[i]]);
  }
  //H298
  fprintf(fptr,"%g\t",JanafEnthalpy(298,&data));
  //S298
  fprintf(fptr,"%g\t",JanafEntropy(298,&data));
  //Cps
  fprintf(fptr,"%g\t",JanafSpecificHeat(300,&data));
  fprintf(fptr,"%g\t",JanafSpecificHeat(400,&data));
  fprintf(fptr,"%g\t",JanafSpecificHeat(500,&data));
  fprintf(fptr,"%g\t",JanafSpecificHeat(600,&data));
  fprintf(fptr,"%g\t",JanafSpecificHeat(700,&data));
  fprintf(fptr,"%g\t",JanafSpecificHeat(800,&data));
  fprintf(fptr,"%g\t",JanafSpecificHeat(1000,&data));
  fprintf(fptr,"%g\n",JanafSpecificHeat(1500,&data));
}

void WriteIsomerJanafBlock(FILE* fptr,
                            int num_print_points,
                            const zerork::ElementalComposition &composition,
                            const std::vector<std::string> &isomer_species_names,
                            const std::vector<JanafThermoData> &isomer_janaf_data)
{
  const double T_epsilon = 1.0e-8;
  double T_min = isomer_janaf_data[0].T_min;
  double T_max = isomer_janaf_data[0].T_max;
  double T_match = isomer_janaf_data[0].T_match;
  double T_plot; //TODO: Fix to plot T_match for all isomers?

  for(size_t i = 1; i < isomer_janaf_data.size(); ++i)
  {
    T_min = std::min(T_min,isomer_janaf_data[i].T_min);
  }
  for(size_t i = 1; i < isomer_janaf_data.size(); ++i)
  {
    T_max = std::max(T_max,isomer_janaf_data[i].T_max);
  }

  fprintf(fptr,"#-----------------------------------------------------------------------------\n");
  fprintf(fptr,"# Composition: %s\n",composition.ToStringWithSeparator("=").c_str());
  fprintf(fptr,"# Species:");
  for(size_t i = 0; i < isomer_species_names.size(); ++i){
    fprintf(fptr,"  %s",isomer_species_names[i].c_str());
  }
  fprintf(fptr,"\n");


  fprintf(fptr, "#   Temp [K]");
  for(size_t i = 0; i < isomer_species_names.size(); ++i) {
    fprintf(fptr,"  Cp/R%12s",isomer_species_names[i].c_str());
  }
  for(size_t i = 0; i < isomer_species_names.size(); ++i){
    fprintf(fptr,"  H/RT%12s",isomer_species_names[i].c_str());
  }
  for(size_t i = 0; i < isomer_species_names.size(); ++i){
    fprintf(fptr,"   S/R%12s",isomer_species_names[i].c_str());
  }
  fprintf(fptr,"\n");

  for(int j=0; j<num_print_points; ++j) {
    T_plot = T_min + 
      (T_match - T_min)*static_cast<double>(j)/(num_print_points-1.0)
      - T_epsilon; // if(T_plot == T_match), then higher branch is used
                   // T_epsilon ensures that a point close to the matching
                   // temperature is used for the lower branch
    fprintf(fptr, "%12.5f",T_plot);
    for(size_t i = 0; i < isomer_species_names.size(); ++i){
        fprintf(fptr, "  %16.7e", JanafSpecificHeatClipped(T_plot,&isomer_janaf_data[i]));
    }
    for(size_t i = 0; i < isomer_species_names.size(); ++i){
        fprintf(fptr, "  %16.7e", JanafEnthalpyClipped(T_plot,&isomer_janaf_data[i]));
    }
    for(size_t i = 0; i < isomer_species_names.size(); ++i){
        fprintf(fptr, "  %16.7e", JanafEntropyClipped(T_plot,&isomer_janaf_data[i]));
    }
    fprintf(fptr,"\n");
  }
  for(int j=0; j<num_print_points; ++j) {
    T_plot = T_match + 
      (T_max - T_match)*static_cast<double>(j)/(num_print_points-1.0);
    fprintf(fptr, "%12.5f",T_plot);
    for(size_t i = 0; i < isomer_species_names.size(); ++i){
        fprintf(fptr, "  %16.7e", JanafSpecificHeatClipped(T_plot,&isomer_janaf_data[i]));
    }
    for(size_t i = 0; i < isomer_species_names.size(); ++i){
        fprintf(fptr, "  %16.7e", JanafEnthalpyClipped(T_plot,&isomer_janaf_data[i]));
    }
    for(size_t i = 0; i < isomer_species_names.size(); ++i){
        fprintf(fptr, "  %16.7e", JanafEntropyClipped(T_plot,&isomer_janaf_data[i]));
    }
    fprintf(fptr,"\n");
  }
  fprintf(fptr,"\n\n");
}
