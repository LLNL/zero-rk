#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros

#ifdef SUNDIALS2
#include <cvode/cvode_dense.h>      // prototypes for CVDense
#include <cvode/cvode_spgmr.h>      // prototypes & constants for CVSPGMR
#elif SUNDIALS3
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#elif SUNDIALS4
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif


#include "user_functions.h"
#include "complete_solver.h"

typedef struct
{
  double fuel_fraction;
  CompositionMap oxidizer;
  double initial_pressure;
  double initial_temperature;
  std::string volume_filename;
  int thist_output;
} InitialConditions;

static int GetNextInitialConditions(const VariableVolumeBatchIFP &parser,
                                    std::ifstream *ic_fptr,
                                    InitialConditions *ic);



int main(int argc, char *argv[])
{
  if(argc != 2) {
    printf("ERROR: incorrect command line usage.\n");
    printf("       use instead %s <input file>\n",argv[0]);
    fflush(stdout);
    exit(-1);
  }
  int error_flag;
  UserData *user_data;
  N_Vector cvode_state;
  std::string thist;
  TimeHistoryParams thist_params;
  void *cvode_memory;
  aux_cvode_structs_t aux_cvode;
  VariableVolumeBatchIFP parser(argv[1]);
  InitialConditions current_ic;
  std::ifstream initial_conditions_file(parser.initialConditionsFile().c_str(),
                                        std::ios_base::in);
  double results[10];
  double simulation_time;
  std::vector<double> initial_mole_fractions;
  std::vector<int> tracked_species_ids;
  FILE *fptr = fopen(parser.outputFile().c_str(),"w");

  if(fptr == NULL) {
    printf("ERROR: can not open output file %s for write operations.\n",
           parser.outputFile().c_str());
    exit(-1);
  }

  // write header
  fprintf(fptr,"# Column  1: [K] initial temperature\n");
  fprintf(fptr,"# Column  2: [Pa] initial pressure\n");
  fprintf(fptr,"# Column  3: [-] mole fraction of composition specified by input file fuelComp map\n");
  fprintf(fptr,"# Column  4: [string] volume time history file\n");
  fprintf(fptr,"# Column  5: [#] CVode error flag (negative values -> error, positive values -> max steps reached)\n");
  fprintf(fptr,"# Column  6: [s] elapsed time from recorded min volume to max dp/dt\n");
  fprintf(fptr,"# Column  7: [s] time at min volume in time history file\n");
  fprintf(fptr,"# Column  8: [m^3] min volume in time history file\n");
  fprintf(fptr,"# Column  9: [s] recorded time at min volume\n");
  fprintf(fptr,"# Column 10: [m^3] recorded min volume\n");
  fprintf(fptr,"# Column 11: [Pa] pressure at recorded min volume\n");
  fprintf(fptr,"# Column 12: [K] temperature at recorded min volume\n");
  fprintf(fptr,"# Column 13: [s] time at max dp/dt\n");
  fprintf(fptr,"# Column 14: [Pa/s] max dp/dt\n");
  fprintf(fptr,"# Column 15: [s] final integrator time\n");
  fprintf(fptr,"# Column 16: [s] CPU wall clock time\n");
  fprintf(fptr,"# Columns 17, ... : Initial mole fractions\n");

  int task_count = 0;
  while(GetNextInitialConditions(parser,
                                 &initial_conditions_file,
                                 &current_ic) != 0) {
    //printf("initial pressure = %.18g [Pa]\n",current_ic.initial_pressure);
    //printf("initial temperature = %.18g [K]\n",current_ic.initial_temperature);
    //printf("volume file = %s\n",current_ic.volume_filename.c_str());

    user_data = new UserData(task_count++,
                             argv[1],
                             current_ic.volume_filename.c_str(),
                             current_ic.initial_pressure,
                             current_ic.initial_temperature,
                             current_ic.fuel_fraction,
                             current_ic.oxidizer);

    if(user_data != NULL) {
      thist_params.record_history    = current_ic.thist_output > 0;
      thist_params.record_prehistory = current_ic.thist_output == 2;
      thist_params.echo_stdout       = false;
      thist_params.step_print_period = user_data->GetParser()->timeHistoryOutputStepPeriod();//10;
      thist_params.min_time_step     = user_data->GetParser()->timeHistoryOutputMinimumTimeStep();//1.0e-5;


      cvode_state = N_VNew_Serial(user_data->GetNumStates());
      user_data->GetInitialState(NV_DATA_S(cvode_state));
#ifdef SUNDIALS4
      cvode_memory = CVodeCreate(CV_BDF);
#else
      cvode_memory = CVodeCreate(CV_BDF, CV_NEWTON);
#endif
      if (CheckFlag((void *)cvode_memory, "CVodeCreate", 0)) {
        exit(-1);
      }

      SetupCompleteSolver(cvode_memory,
                          aux_cvode,
                          user_data,
                          cvode_state,
                          stdout);
      error_flag = SolveVariableVolume(cvode_memory,
                                       user_data,
                                       cvode_state,
                                       NULL,           // A-Factor multipliers
                                       thist_params,   // print out time steps
                                       &thist,
                                       results,
                                       &simulation_time);
      if(error_flag != 0) {
        printf("ERROR: SolveVariableVolume(...) returned error code %d\n",
               error_flag);
        fflush(stdout);
      }
      fprintf(fptr,"%14.7e  %14.7e  %14.7e %s %d  %14.7e",
              current_ic.initial_temperature,
              current_ic.initial_pressure,
              current_ic.fuel_fraction,
              current_ic.volume_filename.c_str(),
              error_flag,
              results[6]-results[2]);
      for(size_t j=0; j<9; ++j) {
        fprintf(fptr, "  %14.7e", results[j]);
      }
      fprintf(fptr,"  %14.7e", simulation_time);
      user_data->GetInitialState(NV_DATA_S(cvode_state));
      user_data->GetTrackedSpeciesIds(tracked_species_ids);
      user_data->GetTrackedMoleFractions(NV_DATA_S(cvode_state),
                                         initial_mole_fractions);
      for(size_t j=0; j<tracked_species_ids.size(); ++j) {
        fprintf(fptr,"  [%s %d] %8.6f",
          user_data->GetReactor()->GetNameOfStateId(tracked_species_ids[j]),
	  tracked_species_ids[j],
	  initial_mole_fractions[j]);
      }

      if(thist_params.record_history) {
        std::string thist_file_name = user_data->GetParser()->timeHistoryOutputFileBase();
        char buffer[48];
        snprintf(buffer,48,"%06d.thist",user_data->GetTaskNum()+1);
        thist_file_name += std::string(buffer);
        FILE *fptr = fopen(thist_file_name.c_str(),"w");
        fprintf(fptr,thist.c_str());
        fclose(fptr);
      }

      fprintf(fptr,"\n");
      fflush(fptr);
      N_VDestroy_Serial(cvode_state);
      CVodeFree(&cvode_memory);
      DestroyAuxCVodeStructs(aux_cvode);
      delete user_data;
    }

  }


  // cvode_state = N_VNew_Serial(user_data->GetNumStates());
  // user_data->GetInitialState(NV_DATA_S(cvode_state));
  // cvode_memory = CVodeCreate(CV_BDF, CV_NEWTON);
  // if (CheckFlag((void *)cvode_memory, "CVodeCreate", 0)) {
  //   exit(-1);
  // }

  // SetupCompleteSolver(cvode_memory,
  //                     user_data,
  //                     cvode_state);

  // double results[3];
  // double simulation_time;
  // printf("# Starting high-level wrapper: SolveVariableVolume\n");
  // fflush(stdout);
  // error_flag = SolveVariableVolume(cvode_memory,
  //                                  user_data,
  //                                  cvode_state,
  //                                  NULL,
  //                                  thist_params,   // print out time steps
  //                                  &results_str,
  //                                  results,
  //                                  &simulation_time);
  // if(error_flag != 0) {
  //   printf("ERROR: SolveVariableVolume(...) returned error code %d\n",
  //          error_flag);
  //   fflush(stdout);
  // }
  // printf("Returned results string:\n %s\n",results_str.c_str());

  // // clean up cvode
  // N_VDestroy_Serial(cvode_state);
  // CVodeFree(&cvode_memory);
  fclose(fptr);
  return 0;
}

template < class ContainerT >
void tokenize(const std::string& str, ContainerT& tokens,
              const std::string& delimiters = " ", bool trimEmpty = false)
{
   if( str.size() == 0 ) return;
   std::string::size_type pos, lastPos = 0, length = str.length();


   while(lastPos < length + 1) {
      pos = str.find_first_of(delimiters, lastPos);
      if(pos == std::string::npos) {
         pos = length;
      }

      if(pos != lastPos || !trimEmpty) {
         tokens.push_back(typename ContainerT::value_type(str.data()+lastPos,
               (typename ContainerT::size_type) pos-lastPos ));
      }

      lastPos = pos + 1;
   }
}


int GetNextInitialConditions(const VariableVolumeBatchIFP &parser,
                             std::ifstream *ic_fptr,
                             InitialConditions *ic)
{
  const char ignore_chars[]="!#";
  size_t ignore_pos;
  std::string line;
  std::string sub_line;

  std::vector<std::string> line_tokens;
  size_t num_oxid_species = parser.oxidizerSpecies().size();

  ic->oxidizer.clear();
  for(size_t i = 0; i < num_oxid_species; ++i) {
    ic->oxidizer[parser.oxidizerSpecies()[i]] = 0.0;
  }

  if(!ic_fptr->is_open()) {
    printf("ERROR: In GetNextInitialConditions(...),\n");
    printf("       initial condition file %s not opened\n",
           parser.initialConditionsFile().c_str());
    return 0;
  }

  if(std::getline(*ic_fptr,line)) {
    ignore_pos = line.find_first_of(ignore_chars);
    sub_line = line.substr(0,ignore_pos);

    line_tokens.clear();
    tokenize<std::vector<std::string> >(sub_line, line_tokens, " \t", true);
    if(line_tokens.size() >= num_oxid_species + 4) {
      ic->fuel_fraction = atof(line_tokens[0].c_str());
      for(size_t i = 0; i < num_oxid_species; ++i) {
        ic->oxidizer[parser.oxidizerSpecies()[i]] = atof(line_tokens[i+1].c_str());
      }
      NormalizeCompositionMap(&(ic->oxidizer));
      ic->initial_temperature = atof(line_tokens[num_oxid_species+1].c_str()) + 273.15;
      ic->initial_pressure = atof(line_tokens[num_oxid_species+2].c_str()) * 1.0e+5;
      std::string filename = line_tokens[num_oxid_species+3];
      ic->volume_filename =
        parser.volumeFilePath() + std::string("/") + filename;

      ic->thist_output = 0;
      if(line_tokens.size() > num_oxid_species + 4) {
        ic->thist_output = atoi(line_tokens[num_oxid_species+4].c_str());
      }

      return 1;
    }
    else {

      return GetNextInitialConditions(parser,
                                      ic_fptr,
                                      ic);
    }
  } else {
    return 0;
  }

}
