#include <stdlib.h>
#include <stdio.h>

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
#endif


#include "user_functions.h"
#include "complete_solver.h"

int main(int argc, char *argv[])
{
  if(argc != 2) {
    printf("ERROR: incorrect command line usage.\n");
    printf("       use instead %s <input file>\n",argv[0]);
    fflush(stdout);
    exit(-1);
  }
  int error_flag;
  UserData user_data(argv[1]);
  N_Vector cvode_state;
  std::string results_str;
  void *cvode_memory;
  aux_cvode_structs_t aux_cvode;
  TimeHistoryParams thist_params;

  // TODO: set the time history params from the file
  thist_params.record_history    = true;
  thist_params.echo_stdout       = true;
  thist_params.step_print_period = user_data.GetParser()->timeHistoryOutputStepPeriod();//10;
  thist_params.min_time_step     = user_data.GetParser()->timeHistoryOutputMinimumTimeStep();//1.0e-5;

  cvode_state = N_VNew_Serial(user_data.GetNumStates());
  user_data.GetInitialState(NV_DATA_S(cvode_state));
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
                      &user_data,
                      cvode_state);

  double results[3];
  double simulation_time;
  printf("# Starting high-level wrapper: SolveVariableVolume\n");
  fflush(stdout);
  error_flag = SolveVariableVolume(cvode_memory,
                                   &user_data,
                                   cvode_state,
                                   NULL,
                                   thist_params,   // print out time steps
                                   &results_str,
                                   results,
                                   &simulation_time);
  if(error_flag != 0) {
    printf("ERROR: SolveVariableVolume(...) returned error code %d\n",
           error_flag);
    fflush(stdout);
  }
  //printf("Returned results string:\n %s\n",results_str.c_str());

  // clean up cvode
  N_VDestroy_Serial(cvode_state);
  CVodeFree(&cvode_memory);
  DestroyAuxCVodeStructs(aux_cvode);

  return 0;
}
