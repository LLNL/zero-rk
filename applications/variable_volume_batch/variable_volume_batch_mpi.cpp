#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "mpi.h"

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
#endif

#include "user_functions.h"
#include "complete_solver.h"

const static int READY_TAG = 1111;
const static int RESULT_TAG = 2222;
const static int TASK_TAG = 3333;

typedef struct
{
  double fuel_fraction;
  CompositionMap oxidizer;
  double initial_pressure;
  double initial_temperature;
  std::string volume_filename;
  int thist_output;
} InitialConditions;

typedef struct
{
  int error_flag;
  double results[10];
  double simulation_time;
  std::vector<double> initial_mole_fractions;
} VariableVolumeResult;

static int GetNextInitialConditions(const VariableVolumeBatchIFP &parser,
                                    std::ifstream *ic_fptr,
                                    InitialConditions *ic);


static int root(int argc, char* argv[]);
static int leaves(int argc, char* argv[]);

static void sendTask(int i, InitialConditions& current_ic, int dest);
static void recieveTask(InitialConditions* current_ic, int* task_num);
static void sendResult(int task_num, int error_flag, std::string& thist, double results[10], double simulation_time, UserData* user_data);
static void recieveResult(VariableVolumeResult* result, int source);
static void printConditions(const InitialConditions& ic);
static void printTaskStatus(int n_completed, int n_total, double elapsed_time);


using zerork::getHighResolutionTime;


int main(int argc, char *argv[])
{
  int rank, status;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  status = 0;
  if(rank==0) {
    if(argc != 2) {
      printf("ERROR: incorrect command line usage.\n");
      printf("       use instead %s <input file>\n",argv[0]);
      fflush(stdout);
      status = 1;
    }
  }
  MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(status != 0) {
    exit(status);
  }

  if(rank == 0) {
    status = root(argc, argv);
  } else {
    status = leaves(argc, argv);
  }

  MPI_Finalize();
  return status;
}

int root(int argc, char* argv[]) {

  int n_procs, n_workers;
  // get worker pool
  MPI_Comm_size(MPI_COMM_WORLD,&n_procs);
  n_workers=n_procs - 1;
  if(n_workers < 1) {
    printf("ERROR: parallel version must have at least 1 worker process.\n");
    printf("       Exiting now.\n");
    exit(-1);
  }

  VariableVolumeBatchIFP parser(argv[1]);
  InitialConditions current_ic;
  std::ifstream initial_conditions_file(parser.initialConditionsFile().c_str(),
                                        std::ios_base::in);
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


  std::vector<InitialConditions> ic_vector;
  while(GetNextInitialConditions(parser,
                                 &initial_conditions_file,
                                 &current_ic) != 0) {
    ic_vector.push_back(current_ic);
  }

  double start_time = getHighResolutionTime();
  int n_in_process = 0;
  int n_tasks = ic_vector.size();
  MPI_Status status;
  std::vector<VariableVolumeResult> result_vector(ic_vector.size());
  for(int i = 0; i < n_tasks; ++i) {
    int statusMsg;
    //Get the reay signal
    MPI_Recv(&statusMsg, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if (status.MPI_TAG == RESULT_TAG) {
      recieveResult(&result_vector[statusMsg],status.MPI_SOURCE);
      n_in_process -= 1;
      printTaskStatus(i-n_in_process,n_tasks,getHighResolutionTime()-start_time);
    }
    sendTask(i, ic_vector[i], status.MPI_SOURCE);
    n_in_process += 1;
  }

  //All tasks submitted.  Get last results and send shutdown signal.
  assert(n_in_process <= n_workers);
  for(int i = 0; i < n_workers; ++i) {
    int statusMsg;
    //Get the reay signal
    MPI_Recv(&statusMsg, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    if (status.MPI_TAG == RESULT_TAG) {
      recieveResult(&result_vector[statusMsg],status.MPI_SOURCE);
      n_in_process -= 1;
      printTaskStatus(n_tasks-n_in_process,n_tasks,getHighResolutionTime()-start_time);
    }
    int done = -1;
    MPI_Send(&done, 1, MPI_INT, status.MPI_SOURCE, TASK_TAG, MPI_COMM_WORLD);
  }
  assert(n_in_process == 0);

  //Process results
  // Get list of species names that we track so we don't have to communicate them;
  std::vector<int> tracked_species_ids;
  std::vector<std::string> tracked_species_names;
  {
    UserData user_data = UserData(-1, argv[1], ic_vector[0].volume_filename.c_str(), ic_vector[0].initial_pressure,
                                           ic_vector[0].initial_temperature, ic_vector[0].fuel_fraction, ic_vector[0].oxidizer);
    user_data.GetTrackedSpeciesIds(tracked_species_ids);
    for(size_t i = 0; i < tracked_species_ids.size(); ++i) {
      tracked_species_names.push_back(user_data.GetReactor()->GetNameOfStateId(tracked_species_ids[i]));
    }
  }

  for(size_t i = 0; i < result_vector.size(); ++i)
  {
    VariableVolumeResult current_result = result_vector[i];

    if(current_result.error_flag != 0) {
      printf("ERROR: SolveVariableVolume(...) returned error code %d\n",
             current_result.error_flag);
      fflush(stdout);
    }
    fprintf(fptr,"%14.7e  %14.7e  %14.7e %s %d  %14.7e",
            ic_vector[i].initial_temperature,
            ic_vector[i].initial_pressure,
            ic_vector[i].fuel_fraction,
            ic_vector[i].volume_filename.c_str(),
            current_result.error_flag,
            current_result.results[6]-current_result.results[2]);
    for(size_t j=0; j<9; ++j) {
      fprintf(fptr, "  %14.7e", current_result.results[j]);
    }
    fprintf(fptr,"  %14.7e", current_result.simulation_time);
    for(size_t j=0; j<tracked_species_ids.size(); ++j) {
      fprintf(fptr,"  [%s %d] %8.6f",
        tracked_species_names[j].c_str(),
        tracked_species_ids[j],
        current_result.initial_mole_fractions[j]);
    }
    fprintf(fptr,"\n");
    fflush(fptr);
  }

  fclose(fptr);
  return 0;
}

int leaves(int argc, char* argv[]) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //We're ready
  int statusMsg = 0;
  MPI_Send(&statusMsg, 1, MPI_INT, 0, READY_TAG, MPI_COMM_WORLD);

  while(1) {
    InitialConditions current_ic;
    //Recieve work
    int task_num;
    recieveTask(&current_ic, &task_num);
    if(task_num < 0) break;

    UserData* user_data = new UserData(task_num,
                              argv[1],
                              current_ic.volume_filename.c_str(),
                              current_ic.initial_pressure,
                              current_ic.initial_temperature,
                              current_ic.fuel_fraction,
                              current_ic.oxidizer);

    if(user_data != NULL) {
      int error_flag;
      N_Vector cvode_state;
      std::string thist;
      void *cvode_memory;
      aux_cvode_structs_t aux_cvode;
      double results[10];
      double simulation_time;
      TimeHistoryParams thist_params;

      // TODO: set the time history params from the file
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

      if(thist_params.record_history) {
        std::string thist_file_name = user_data->GetParser()->timeHistoryOutputFileBase();
        char buffer[48];
        snprintf(buffer,48,"%06d.thist",user_data->GetTaskNum()+1);
        thist_file_name += std::string(buffer);
        FILE *fptr = fopen(thist_file_name.c_str(),"w");
        fprintf(fptr,thist.c_str());
        fclose(fptr);
      }

      sendResult(task_num, error_flag, thist, results, simulation_time, user_data);
      N_VDestroy_Serial(cvode_state);
      CVodeFree(&cvode_memory);
      DestroyAuxCVodeStructs(aux_cvode);

      delete user_data;
    }
  }

  return 0;
}

void sendTask(int i, InitialConditions& current_ic, int dest) {
  MPI_Send(&i, 1, MPI_INT, dest, TASK_TAG, MPI_COMM_WORLD);
  MPI_Send(&current_ic.fuel_fraction, 1, MPI_DOUBLE, dest, TASK_TAG, MPI_COMM_WORLD);

  CompositionMap& oxidizer = current_ic.oxidizer;
  int n_oxidizer_species = oxidizer.size();
  MPI_Send(&n_oxidizer_species, 1, MPI_INT, dest, TASK_TAG, MPI_COMM_WORLD);
  for(CompositionMap::const_iterator it=oxidizer.begin();
      it!=oxidizer.end(); ++it)
  {
    std::string species = it->first;
    double value = it->second;
    MPI_Send(&(species[0]), species.size(), MPI_CHAR, dest, TASK_TAG, MPI_COMM_WORLD);
    MPI_Send(&value, 1, MPI_DOUBLE, dest, TASK_TAG, MPI_COMM_WORLD);
  }

  MPI_Send(&current_ic.initial_pressure, 1, MPI_DOUBLE, dest, TASK_TAG, MPI_COMM_WORLD);
  MPI_Send(&current_ic.initial_temperature, 1, MPI_DOUBLE, dest, TASK_TAG, MPI_COMM_WORLD);
  MPI_Send(&(current_ic.volume_filename[0]), current_ic.volume_filename.size(), MPI_CHAR, dest, TASK_TAG, MPI_COMM_WORLD);
  MPI_Send(&current_ic.thist_output, 1, MPI_INT, dest, TASK_TAG, MPI_COMM_WORLD);
}

void recieveTask(InitialConditions* current_ic, int* task_num) {
  MPI_Status status;
  MPI_Recv(task_num, 1, MPI_INT, 0, TASK_TAG, MPI_COMM_WORLD, &status);
  if(*task_num < 0) return;
  MPI_Recv(&(current_ic->fuel_fraction), 1, MPI_DOUBLE, 0, TASK_TAG, MPI_COMM_WORLD, &status);

  current_ic->oxidizer.clear();
  int num_oxidizer_species;
  MPI_Recv(&num_oxidizer_species, 1, MPI_INT, 0, TASK_TAG, MPI_COMM_WORLD, &status);
  for(int i = 0; i < num_oxidizer_species; ++i)
  {
    MPI_Probe(0, TASK_TAG, MPI_COMM_WORLD, &status);
    int num_chars;
    MPI_Get_count(&status, MPI_CHAR, &num_chars);

    double value;
    std::string species(" ", num_chars);
    MPI_Recv(&(species[0]), num_chars, MPI_CHAR, 0, TASK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&value, 1, MPI_DOUBLE, 0, TASK_TAG, MPI_COMM_WORLD, &status);
    current_ic->oxidizer[species] = value;
  }

  MPI_Recv(&(current_ic->initial_pressure), 1, MPI_DOUBLE, 0, TASK_TAG, MPI_COMM_WORLD, &status);
  MPI_Recv(&(current_ic->initial_temperature), 1, MPI_DOUBLE, 0, TASK_TAG, MPI_COMM_WORLD, &status);

  MPI_Probe(0, TASK_TAG, MPI_COMM_WORLD, &status);
  int num_chars;
  MPI_Get_count(&status, MPI_CHAR, &num_chars);

  current_ic->volume_filename.resize(num_chars);
  MPI_Recv(&(current_ic->volume_filename[0]), num_chars, MPI_CHAR, 0, TASK_TAG, MPI_COMM_WORLD, &status);

  MPI_Recv(&(current_ic->thist_output), 1, MPI_INT, 0, TASK_TAG, MPI_COMM_WORLD, &status);
}


void sendResult(int task_num, int error_flag, std::string& thist, double results[10], double simulation_time, UserData* user_data)
{
  std::vector<double> initial_mole_fractions;
  std::vector<int> tracked_species_ids;
  N_Vector cvode_state = N_VNew_Serial(user_data->GetNumStates());
  user_data->GetInitialState(NV_DATA_S(cvode_state));
  user_data->GetTrackedSpeciesIds(tracked_species_ids);
  user_data->GetTrackedMoleFractions(NV_DATA_S(cvode_state),
                                     initial_mole_fractions);

  N_VDestroy_Serial(cvode_state);
  MPI_Send(&task_num, 1, MPI_INT, 0, RESULT_TAG, MPI_COMM_WORLD);

  MPI_Send(&error_flag, 1, MPI_INT, 0, RESULT_TAG, MPI_COMM_WORLD);
  MPI_Send(results, 10, MPI_DOUBLE, 0, RESULT_TAG, MPI_COMM_WORLD);
  MPI_Send(&simulation_time, 1, MPI_DOUBLE, 0, RESULT_TAG, MPI_COMM_WORLD);

  int n_mole_fractions = initial_mole_fractions.size();
  MPI_Send(&n_mole_fractions, 1, MPI_INT, 0, RESULT_TAG, MPI_COMM_WORLD);
  MPI_Send(&(initial_mole_fractions[0]), n_mole_fractions, MPI_DOUBLE, 0, RESULT_TAG, MPI_COMM_WORLD);
}

void recieveResult(VariableVolumeResult* result, int source)
{
  MPI_Status status;
  MPI_Recv(&(result->error_flag), 1, MPI_INT, source, RESULT_TAG, MPI_COMM_WORLD, &status);
  MPI_Recv(result->results, 10, MPI_DOUBLE, source, RESULT_TAG, MPI_COMM_WORLD, &status);
  MPI_Recv(&(result->simulation_time), 1, MPI_DOUBLE, source, RESULT_TAG, MPI_COMM_WORLD, &status);
  int n_mole_fractions;
  MPI_Recv(&n_mole_fractions, 1, MPI_INT, source, RESULT_TAG, MPI_COMM_WORLD, &status);
  result->initial_mole_fractions.resize(n_mole_fractions);
  MPI_Recv(&(result->initial_mole_fractions[0]), n_mole_fractions, MPI_DOUBLE, source, RESULT_TAG, MPI_COMM_WORLD, &status);
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

void printConditions(const InitialConditions& ic) {
  printf("fuel_fraction: %20.13g\n",ic.fuel_fraction);
  printf("oxodizer:");
  for(CompositionMap::const_iterator it=ic.oxidizer.begin();
      it!=ic.oxidizer.end(); ++it)
  {
    std::string species = it->first;
    double value = it->second;
    printf(" [%s] %20.13g",species.c_str(), value);
  }
  printf("\ninitial_pressure: %20.13g\n",ic.initial_pressure);
  printf("initial_temperature: %20.13g\n",ic.initial_temperature);
  printf("volume_filename: %s\n",ic.volume_filename.c_str());
}

void printTaskStatus(int n_completed, int n_total, double elapsed_time) {
  if(n_completed > 0) {
    double avg_time_per_task = elapsed_time/n_completed;
    double est_time_left = (n_total-n_completed)*avg_time_per_task;
    printf("Task status: %d/%d  completed in %8.3g seconds.  Estimated time remaining: %8.3g seconds\n",
           n_completed, n_total, elapsed_time, est_time_left);
    fflush(stdout);
  }
}
