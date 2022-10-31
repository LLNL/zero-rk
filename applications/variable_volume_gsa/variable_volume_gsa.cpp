#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

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

#include <mechanism_info/mechanism_info.h>

#include "utilities/file_utilities.h"
#include "utilities/sequential_file_matrix.h"
#include "gsa_stats.h"

#include "user_functions_gsa.h"
#include "complete_solver_gsa.h"

using zerork::getHighResolutionTime;

const int MASTER_RANK=0;
const int WORK_TAG=1;
const int KILL_TAG=2;
const double SHUT_DOWN_TIME = 2.0;

int MasterTasks(int inp_argc, char *inp_argv[]);
int WorkerTasks(int inp_argc, char *inp_argv[]);
void PrintEstimateEndTime(const int num_done,
                          const int num_total,
                          const double start_time);
void SendWorkersKillTag();


class FlushableIndexList
{
 public:
  FlushableIndexList(const size_t list_size);
  ~FlushableIndexList() {};
  size_t Flush(std::vector<size_t> *index_list);
  // Flush returns the size of the index list that was flushed
  bool Ready(const size_t index);
  // returns true if the index was set to the ready for flush state, it will
  // return false if the index was out-of-bounds

 private:
  size_t list_size_;
  size_t next_index_;
  std::vector<bool> ready_;
  std::vector<bool> flush_;
};


int main(int argc, char *argv[])
{
  int my_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if(my_rank == MASTER_RANK) {
    MasterTasks(argc, argv);
  }
  else {
    WorkerTasks(argc, argv);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
int MasterTasks(int inp_argc, char *inp_argv[])
{
  double start_time = getHighResolutionTime();
  double elapsed_time, sum_worker_time;
  int num_threads;
  int num_workers;
  int num_tasks;

  int task_id;
  int task_id_tag;
  int recv_count=0;

  MPI_Status status;

  // worker_solution[0] = final ODE time
  // worker_solution[1] = max dP/dt
  // worker_solution[2] = time of max dP/dt
  // worker_solution[3] = total chemical heat release
  // worker_solution[msg_send_size-1] = CPU wall clock time
  int msg_send_size;
  int num_burn_fractions;
  int num_hrr_powers;
  double *worker_solution;
  double *master_solution;
  double *all_solutions;

  int num_reactions=0;
  double *afactor_mult;

  int error_flag;
  N_Vector cvode_state;
  std::string results_str;
  void *cvode_memory=NULL;
  aux_cvode_structs_t aux_cvode;
  UserData *user_data=NULL;
  zerork::utilities::SequentialFileMatrixRead *gsa_matrix=NULL;
  GSAStats *stats=NULL;
  TimeHistoryParams thist_params;

  double ref_heat_release;

  FILE *thist_fptr;
  FILE *idt_output_fptr;
  FILE *gsa_check_fptr;

  if(inp_argc != 2) {
    printf("ERROR: Incorrect command line usage,\n");
    printf("       use instead %s <input_file_name>\n",inp_argv[0]);
    printf("       Exiting application in %6.2f seconds.\n",SHUT_DOWN_TIME);
    fflush(stdout);
    SendWorkersKillTag();
    return -1;
  }
  if(!zerork::utilities::FileIsReadable(inp_argv[1])) {
    printf("ERROR: Could not open input file %s for read operation\n",
           inp_argv[1]);
    printf("       Exiting application in %6.2f seconds.\n",SHUT_DOWN_TIME);
    fflush(stdout);
    SendWorkersKillTag();
    return -2;
  }
  user_data = new UserData(inp_argv[1], // input filename
                           true);      // don't write mechanism parser log
  num_burn_fractions =
    static_cast<int>(user_data->GetParser()->burnFraction().size());
  num_hrr_powers =
    static_cast<int>(user_data->GetParser()->hrrIntegralPower().size());

  // add one to the results size for the CPU time
  msg_send_size = NUM_FIXED_RESULTS + 1 + num_burn_fractions + num_hrr_powers;

  worker_solution = new double[msg_send_size];
  master_solution = new double[msg_send_size];
  //printf("msg_send_size = %d\n",msg_send_size);


  if(user_data == NULL) {
    printf("ERROR: UserData constructor failed for input file %s.\n",
           inp_argv[1]);
    printf("       Exiting application in %6.2f seconds.\n",SHUT_DOWN_TIME);
    fflush(stdout);
    SendWorkersKillTag();
    return -3;
  }
  num_reactions = user_data->GetReactor()->GetNumReactions();
  afactor_mult = new double[num_reactions];

  gsa_matrix =
    new zerork::utilities::SequentialFileMatrixRead(user_data->GetParser()->gsaMatrixFile().c_str());
  if(gsa_matrix == NULL) {
    printf("ERROR: SequentialFileMatrixRead constructor failed for input file %s.\n",
           user_data->GetParser()->gsaMatrixFile().c_str());
    printf("       Exiting application in %6.2f seconds.\n",SHUT_DOWN_TIME);
    fflush(stdout);
    SendWorkersKillTag();
    return -4;
  }
  if(gsa_matrix->num_columns() != num_reactions) {
    printf("ERROR: GSA Matrix file is incompatible with the mechanism.\n");
    printf("       The number of matrix columns must match the number of reactions.\n");
    printf("         GSA Matrix file %s number of columns   = %d\n",
           gsa_matrix->filename(),gsa_matrix->num_columns());
    printf("         Mechanism file %s  number of reactions = %d\n",
           user_data->GetParser()->mechFile().c_str(),num_reactions);
    printf("       Exiting application in %6.2f seconds.\n",SHUT_DOWN_TIME);
    fflush(stdout);
    SendWorkersKillTag();
    return -5;
  }
  num_tasks = gsa_matrix->num_rows();
  all_solutions = new double[num_tasks*msg_send_size];
  for(int j=0; j<(num_tasks*msg_send_size); ++j) {
     all_solutions[j] = 0.0;
  }

  FlushableIndexList complete_tasks(num_tasks);
  size_t num_writes;
  std::vector<size_t> write_list;

  stats = new GSAStats(num_reactions);
  if(stats == NULL) {
    printf("ERROR: GSAStats constructor failed for input file %d reactions.\n",
           num_reactions);
    printf("       Exiting application in %6.2f seconds.\n",SHUT_DOWN_TIME);
    fflush(stdout);
    SendWorkersKillTag();
    return -6;
  }

  // initialize state vector for cvode
  cvode_state = N_VNew_Serial(user_data->GetNumStates());

  // set the initial state
  user_data->GetInitialState(NV_DATA_S(cvode_state));
  // set the reference heat release to a nominal value to prevent
  // accessing an undefined value or dividing by zero, the master thread
  // does not use the heat release fraction results
  user_data->SetReferenceHeatRelease(user_data->GetParser()->refHeatRelease());

  // allocate initial memory
#ifdef SUNDIALS4
  cvode_memory = CVodeCreate(CV_BDF);
#else
  cvode_memory = CVodeCreate(CV_BDF, CV_NEWTON);
#endif
  if (CheckFlag((void *)cvode_memory, "CVodeCreate", 0)) {
    printf("ERROR: In MasterTask() (rank=%d),\n",MASTER_RANK);
    printf("       CVodeCreate(CV_BDF, CV_NEWTON) failed.\n");
    printf("       Exiting application in %6.2f seconds.\n",SHUT_DOWN_TIME);
    fflush(stdout);
    SendWorkersKillTag();
    return -7;
  }
  // complete the solver setup
  SetupCompleteSolver(cvode_memory,
                      aux_cvode,
                      user_data,
                      cvode_state);

  MPI_Comm_size(MPI_COMM_WORLD,&num_threads);
  num_workers = num_threads-1;
  if(num_workers < 1) {
    printf("ERROR: In MasterTask() (rank=%d),\n",MASTER_RANK);
    printf("       current version does not support MPI with %d worker threads.\n",num_workers);
    fflush(stdout);
    return -8;
  }
  printf("# Starting unperturbed simulation on master rank...\n");
  fflush(stdout);
  for(int j=0; j<num_reactions; ++j) {
    afactor_mult[j] = 1.0;
  }
  // TODO: set time history record from input file
  thist_params.record_history    = true;
  thist_params.echo_stdout       = user_data->GetParser()->thistEcho();
  thist_params.step_print_period = user_data->GetParser()->thistStepPeriod();
  thist_params.min_time_step     = user_data->GetParser()->thistMinDt();

  error_flag = SolveVariableVolume(cvode_memory,
                                   user_data,
                                   cvode_state,
                                   &afactor_mult[0],
                                   thist_params,   // print out time steps
                                   &results_str,
                                   &master_solution[0],
                                   &master_solution[msg_send_size-1]);
  if(error_flag != 0) {
    printf("ERROR: SolveVariableVolume(...) returned error code %d\n",
           error_flag);
    printf("       for the unperturbed A-Factor case on the master thread.\n");
    printf("       Resolve the errors before running the GSA.\n");
    printf("   Worker rank id             [#]: %d\n",MASTER_RANK);
    printf("   Task id                    [#]: %d\n",0);
    printf("   Final ODE time             [s]: %.18g\n",
           master_solution[0]);
    printf("   Max dP/dt               [Pa/s]: %.18g\n",
           master_solution[1]);
    printf("   Time at max dP/dt          [s]: %.18g\n",
           master_solution[2]);
    printf("   Simulation wall clock time [s]: %.18g\n",
           master_solution[msg_send_size-1]);
    printf("   Exiting application in %6.2f seconds.\n",SHUT_DOWN_TIME);
    fflush(stdout);
    SendWorkersKillTag();
    return -9;
  }
  // set reference heat release
  ref_heat_release = master_solution[3];
  if(user_data->GetParser()->useRefHeatRelease()) {
    ref_heat_release = user_data->GetParser()->refHeatRelease();
  }
  // set the reference heat release from the master thread
  MPI_Bcast(&ref_heat_release,
            1,
            MPI_DOUBLE,
            MASTER_RANK,
            MPI_COMM_WORLD);

  printf("#   Unperturbed A-Factor time of max dP/dt       [s]: %14.7e\n",
         master_solution[2]);
  printf("#   Unperturbed A-Factor total chem heat release [J]: %14.7e\n",
         master_solution[3]);
  printf("#   Reference total chem heat release            [J]: %14.7e\n",
         ref_heat_release);
  printf("#   Single case simulation time                  [s]: %14.7e\n",
         master_solution[msg_send_size-1]);
  printf("#   Estimated remaining wall clock time          [s]: %14.7e\n",
         master_solution[msg_send_size-1]*
         (double)num_tasks/(double)num_workers);
  printf("# -----------------------------------------------------------------------------\n");
  fflush(stdout);
  sum_worker_time=master_solution[msg_send_size-1];

  // clean up master cvode information
  N_VDestroy_Serial(cvode_state);
  CVodeFree(&cvode_memory);
  DestroyAuxCVodeStructs(aux_cvode);
  // -------------------------------------------------------------------------
  // record the results from the unperturbed solution
  // -------------------------------------------------------------------------
  // setup header for the results collected
  idt_output_fptr = fopen(user_data->GetParser()->outFile().c_str(),"w");
  if (idt_output_fptr == NULL) {
    printf("ERROR: In MasterTask() (rank=%d),\n",MASTER_RANK);
    printf("       Failed to create output file: %s.\n", user_data->GetParser()->outFile().c_str());
    printf("       Exiting application in %6.2f seconds.\n",SHUT_DOWN_TIME);
    fflush(stdout);
    SendWorkersKillTag();
    return -10;
  }

  fprintf(idt_output_fptr,
          "# Column  1: [#]    task id (zero is the unperturbed solution)\n");
  fprintf(idt_output_fptr,
          "# Column  2: [s]    final time reached by the ODE integrator\n");
  fprintf(idt_output_fptr,
          "# Column  3: [Pa/s] max dP/dt pressure rise rate (1st order backward difference\n");
  fprintf(idt_output_fptr,
          "# Column  4: [s]    time at max dP/dt\n");
  fprintf(idt_output_fptr,
          "# Column  5: [J]    total chemical heat release\n");
  fprintf(idt_output_fptr,
          "# Column  6: [J/s]  max chemical heat release rate (linear)\n");
  fprintf(idt_output_fptr,
          "# Column  7: [J/s]  max chemical heat release rate (quadratic)\n");

  for(int j=0; j<num_hrr_powers; ++j) {
    fprintf(idt_output_fptr,
            "# Column %2d: [J/s]    [integral(heat release rate)^%f ]^(1/%f)\n",
            NUM_FIXED_RESULTS+2+j,
            user_data->heat_release_rate_powers_[j+1],
            user_data->heat_release_rate_powers_[j+1]);
  }
  for(int j=0; j<num_burn_fractions; ++j) {

    fprintf(idt_output_fptr,
            "# Column %d: [s]    time at burn fraction %4.3f\n",
            num_hrr_powers+NUM_FIXED_RESULTS+1+j,user_data->GetParser()->burnFraction()[j]);
  }
  fprintf(idt_output_fptr,
          "# Column %d: [s]    solver wall clock time\n",
	  num_hrr_powers+NUM_FIXED_RESULTS+1+num_burn_fractions);

  fprintf(idt_output_fptr, "%6d", 0);
  for(int j=0; j<(NUM_FIXED_RESULTS+num_hrr_powers); ++j) {
    fprintf(idt_output_fptr, "  %14.7e",master_solution[j]);
  }
  for(int j=0; j<num_burn_fractions; ++j) {
    fprintf(idt_output_fptr,"                ");
    // TODO: print out the burn fractions if the user specified a reference
    // heat release
  }
  fprintf(idt_output_fptr,"  %14.7e\n",master_solution[msg_send_size-1]);
  fflush(idt_output_fptr);
  // assign the first batch of IDT problems to the workers

  if(num_workers > num_tasks) {
    printf("WARNING: the number of worker threads (%d) exceeds\n",
           num_workers);
    printf("         the number of AFactor perturbation tasks (%d).\n",
           num_tasks);
    num_workers=num_tasks;
  }
  task_id = 0;
  // CAUTION: assumes master rank is zero
  for(int j=1; j<=num_workers; j++) {

    gsa_matrix->GetNextRow(&afactor_mult[0]);
    stats->AddNextSample(&afactor_mult[0]);

    MPI_Send(&task_id,         // message buffer (task index)
             1,                // message buffer count
             MPI_INT,          // message buffer type
             j,                // process rank
             WORK_TAG,         // user-specified message tag
             MPI_COMM_WORLD);  // default communicator for all threads

    MPI_Send(&afactor_mult[0], // message buffer
             num_reactions,    // message buffer count
             MPI_DOUBLE,       // message buffer type
             j,                // process rank
             WORK_TAG,         // user-specified message tag
             MPI_COMM_WORLD);  // default communicator for all threads


    ++task_id;
   }
  // -------------------------------------------------------------------------
  // start output files
  thist_fptr = fopen(user_data->GetParser()->thistFile().c_str(),"w");
  fprintf(thist_fptr,"%s",results_str.c_str());
  fclose(thist_fptr);
  gsa_check_fptr  = fopen(user_data->GetParser()->checkFile().c_str(),"w");
  // -------------------------------------------------------------------------
  // receive the results and send out the remaining tasks
  while(task_id < num_tasks) {
    // recieve result from any worker
    MPI_Recv(&worker_solution[0], // results buffer
             msg_send_size,       // size of buffer
             MPI_DOUBLE,          // buffer type
             MPI_ANY_SOURCE,      // receive msgResults from any thread
             MPI_ANY_TAG,         // receive msgResults with any tags
             MPI_COMM_WORLD,      // default communicator for all threads
             &status);            // info about the communication

    ++recv_count;
    if(recv_count%user_data->GetParser()->taskProgressPeriod() == 0) {
      PrintEstimateEndTime(recv_count,num_tasks,start_time);
    }

    // record the results for the task id
    task_id_tag = status.MPI_TAG; // workers send the taskId in the MPI_TAG

    for(int j=0; j<msg_send_size; j++) {
      all_solutions[task_id_tag*msg_send_size+j]=worker_solution[j];
    }
    sum_worker_time+=worker_solution[msg_send_size-1];
    // --------------------------------------------------------------------
    // record that the task is completed and check the number of consecutive
    // completed task indexes that can be flushed
    complete_tasks.Ready(task_id_tag);
    num_writes = complete_tasks.Flush(&write_list);

    for(size_t j=0; j<num_writes; ++j) {

      fprintf(idt_output_fptr, "%6d", (int)write_list[j]+1);
      for(size_t k=0; k<(size_t)msg_send_size; ++k) {
        fprintf(idt_output_fptr,
                "  %14.7e",
                all_solutions[write_list[j]*msg_send_size + k]);
      }
      fprintf(idt_output_fptr,"\n");
    }
    fflush(idt_output_fptr);

    // record and write complete tasks
    // --------------------------------------------------------------------

    // get the next set of A-Factor pertubation multipliers
    gsa_matrix->GetNextRow(&afactor_mult[0]);
    stats->AddNextSample(&afactor_mult[0]);

    // assign the worker a new problem
    MPI_Send(&task_id,          // message buffer (rxn idx to perturb)
             1,                 // message buffer count
             MPI_INT,           // message buffer type
             status.MPI_SOURCE, // process rank that just finished
             WORK_TAG,          // user-specified message tag
             MPI_COMM_WORLD);   // default communicator for all threads


    MPI_Send(&afactor_mult[0],  // message buffer
             num_reactions,     // message buffer count
             MPI_DOUBLE,        // message buffer type
             status.MPI_SOURCE, // process rank that just finished
             WORK_TAG,          // user-specified message tag
             MPI_COMM_WORLD);   // default communicator for all threads

    ++task_id;
  }
  // -------------------------------------------------------------------------
  // All tasks have been assigned, and all workers (rank 1,...,nWorker)
  // received at least one task.  Now collect the last results from each
  // worker.
  for(int j=1; j<=num_workers; j++) {
    MPI_Recv(&worker_solution[0], // results buffer
             msg_send_size,       // size of buffer
             MPI_DOUBLE,          // buffer type
             MPI_ANY_SOURCE,      // receive msgResults from any thread
             MPI_ANY_TAG,         // receive msgResults with any tags
             MPI_COMM_WORLD,      // default communicator for all threads
             &status);            // info about the communication
    ++recv_count;
    if(recv_count%user_data->GetParser()->taskProgressPeriod() == 0) {
      PrintEstimateEndTime(recv_count,num_tasks,start_time);
    }

    // record the results for the task id
    task_id_tag = status.MPI_TAG; // workers send the taskId in the MPI_TAG

    for(int k=0; k<msg_send_size; k++) {
      all_solutions[task_id_tag*msg_send_size+k]=worker_solution[k];
    }
    sum_worker_time+=worker_solution[msg_send_size-1];
    // --------------------------------------------------------------------
    // record that the task is completed and check the number of consecutive
    // completed task indexes that can be flushed
    complete_tasks.Ready(task_id_tag);
    num_writes = complete_tasks.Flush(&write_list);

    for(size_t k=0; k<num_writes; ++k) {

      fprintf(idt_output_fptr, "%6d", (int)write_list[k]+1);
      for(size_t m=0; m<(size_t)msg_send_size; ++m) {
        fprintf(idt_output_fptr,
                "  %14.7e",
                all_solutions[write_list[k]*msg_send_size + m]);
      }
      fprintf(idt_output_fptr,"\n");
    }
    fflush(idt_output_fptr);

    // record and write complete tasks
    // --------------------------------------------------------------------
  }

  // -------------------------------------------------------------------------
  // no more work - send out the worker KILL_TAG with an empty message
  SendWorkersKillTag();

  // Double check there are no remaining completed tasks to flush
  // --------------------------------------------------------------------
  // record that the task is completed and check the number of consecutive
  // completed task indexes that can be flushed
  num_writes = complete_tasks.Flush(&write_list);

  for(size_t k=0; k<num_writes; ++k) {

    fprintf(idt_output_fptr, "%6d", (int)write_list[k]+1);
    for(size_t m=0; m<(size_t)msg_send_size; ++m) {
      fprintf(idt_output_fptr,
              "  %14.7e",
              all_solutions[write_list[k]*msg_send_size + m]);
    }
    fprintf(idt_output_fptr,"\n");
  }
  // record and write complete tasks
  // --------------------------------------------------------------------



  fclose(idt_output_fptr);

  // reaction check for bounds min, max A-Factor and estimate
  // of normal parameters for log(afactor_mult)
  MechanismInfo mech_info(user_data->GetParser()->mechFile().c_str(),
			  user_data->GetParser()->thermFile().c_str(),
                          "");
  fprintf(gsa_check_fptr,
          "# Basic distribution statistics for each reaction's perturbation multiplier.\n");
  fprintf(gsa_check_fptr,
          "# Here 'a' is the perturbation multiplier, so a=1 is unperturbed.\n");
  fprintf(gsa_check_fptr,"# Note reaction id is index-0.\n");
  fprintf(gsa_check_fptr,"# Number of samples   = %d (GSA matrix rows)\n",
          stats->num_samples());
  fprintf(gsa_check_fptr,"# Number of reactions = %d (GSA matrix columns)\n",
          stats->num_dimensions());
  fprintf(gsa_check_fptr,"# -----------------------------------------------------------------------------\n");
  fprintf(gsa_check_fptr,"#   id  reaction name                            min(a)         max(a)   mean[log(a)]  stdev[log(a)]\n");
  for(int j=0; j<num_reactions; ++j) {
    // note that '-' indicates left-justify for value
    fprintf(gsa_check_fptr,"%6d  %-32s  %13.4e  %13.4e  %13.4e  %13.4e\n",
            j,
            //idtCtrl.mech->getReactionName(j),
            mech_info.GetReactionName(j),
            stats->MinValue(j),
            stats->MaxValue(j),
            stats->LogMean(j),
            sqrt(stats->LogVariance(j)));
  }


  elapsed_time = getHighResolutionTime() - start_time;
  printf("# Elapsed wall clock time  [s]: %14.7e (%d worker threads)\n",
         elapsed_time,num_workers);
  printf("# Total CPU time of solver [s]: %14.7e\n",sum_worker_time);
  fflush(stdout);
  // clean up memory on Master

  fclose(gsa_check_fptr);

  delete[] master_solution;
  delete[] worker_solution;

  delete user_data;
  delete gsa_matrix;
  delete stats;
  delete[] afactor_mult;
  delete[] all_solutions;

  return 0;
}

int WorkerTasks(int inp_argc, char *inp_argv[])
{
  int my_rank;
  int my_task_id;
  int num_runs=0;
  double *afactor_mult=NULL;
  MPI_Status status;
  int num_reactions=0;

  // worker_solution[0] = final ODE time
  // worker_solution[1] = max dP/dt
  // worker_solution[2] = time of max dP/dt
  // worker_solution[3] = total heat release
  // worker_solution[msg_send_size-1] = CPU wall clock time
  int msg_send_size;
  double *worker_solution;

  int error_flag;
  N_Vector cvode_state;
  std::string results_str;
  void *cvode_memory;
  aux_cvode_structs_t aux_cvode;
  UserData *user_data=NULL;
  TimeHistoryParams thist_params;

  double ref_heat_release;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // receive the reference heat release from the master thread
  MPI_Bcast(&ref_heat_release,
            1,
            MPI_DOUBLE,
            MASTER_RANK,
            MPI_COMM_WORLD);
  //printf("# worker %d set reference heat release to %14.7e\n",
  //       my_rank,ref_heat_release);

  while(1) {
    MPI_Recv(&my_task_id,     // buffer - task id
             1,               // buffer size
             MPI_INT,         // buffer type
             MASTER_RANK,     // only receive messages from the MASTER
             MPI_ANY_TAG,     // receive messages with any tag
             MPI_COMM_WORLD,  // use the default communicator with all threads
             &status);        // info about the message

    if(status.MPI_TAG == KILL_TAG) {
      //printf("Worker id %d received KILL_TAG\n",my_rank);
      //fflush(stdout);
      break; // leave the loop when the master thread send the KILL_TAG
    }

    if(num_runs == 0) {
      // perform the inital setup
      user_data = new UserData(inp_argv[1], // input filename
                               false);      // don't write mechanism parser log

      user_data->SetReferenceHeatRelease(ref_heat_release);

      num_reactions = user_data->GetReactor()->GetNumReactions();
      afactor_mult = new double[num_reactions];

      // add one the result
      msg_send_size =
        NUM_FIXED_RESULTS+1+static_cast<int>(user_data->GetParser()->burnFraction().size())+static_cast<int>(user_data->GetParser()->hrrIntegralPower().size());
      worker_solution = new double[msg_send_size];

      // initialize state vector for cvode
      cvode_state = N_VNew_Serial(user_data->GetNumStates());

      // set the initial state
      user_data->GetInitialState(NV_DATA_S(cvode_state));

      // allocate initial memory
#ifdef SUNDIALS4
      cvode_memory = CVodeCreate(CV_BDF);
#else
      cvode_memory = CVodeCreate(CV_BDF, CV_NEWTON);
#endif
      if (CheckFlag((void *)cvode_memory, "CVodeCreate", 0)) {
        printf("ERROR: In WorkerTask() (rank=%d),\n",my_rank);
        printf("       CVodeCreate(CV_BDF, CV_NEWTON) failed.\n");
        fflush(stdout);
        exit(-1);
      }
      // complete the solver setup
      SetupCompleteSolver(cvode_memory,
                          aux_cvode,
                          user_data,
                          cvode_state);

      thist_params.record_history    = false;
      thist_params.echo_stdout       = false;
      thist_params.step_print_period = 1000000000;
      thist_params.min_time_step     = 1.0e300;

      thist_params.record_history    = user_data->GetParser()->perturbationTimeHistory();
      if(thist_params.record_history) {
        // Don't echo
        //thist_params.echo_stdout       = user_data->GetParser()->thistEcho();
        thist_params.step_print_period = user_data->GetParser()->thistStepPeriod();
        thist_params.min_time_step     = user_data->GetParser()->thistMinDt();
      }
    } // end the initial setup when (num_runs == 0)

    MPI_Recv(&afactor_mult[0], // message buffer
             num_reactions,     // message buffer count
             MPI_DOUBLE,       // message buffer type
             MASTER_RANK,      // only receive messages from the MASTER
             WORK_TAG,         // user-specified message tag
             MPI_COMM_WORLD,   // default communicator for all threads
             &status);

    // solve the ignition delay for the conditions
    error_flag = SolveVariableVolume(cvode_memory,
                                     user_data,
                                     cvode_state,
                                     &afactor_mult[0],
                                     thist_params,   // print out time steps
                                     &results_str,
                                     &worker_solution[0],
                                     &worker_solution[msg_send_size-1]);

    ++num_runs;
    if(error_flag != 0) {
      printf("WARNING: SolveVariableVolume(...) returned error code %d\n",
             error_flag);
      printf("   Worker rank id             [#]: %d\n",my_rank);
      printf("   Task id                    [#]: %d\n",my_task_id);
      printf("   Worker solve count         [#]: %d\n",num_runs);
      printf("   Final ODE time             [s]: %.18g\n",
             worker_solution[0]);
      printf("   Max dP/dt               [Pa/s]: %.18g\n",
             worker_solution[1]);
      printf("   Time at max dP/dt          [s]: %.18g\n",
             worker_solution[2]);
      printf("   Simulation wall clock time [s]: %.18g\n",
             worker_solution[msg_send_size-1]);

      fflush(stdout);
    }
    if(thist_params.record_history) {
      char my_thist_file_name[255];
      snprintf(my_thist_file_name,255,"%s_%08d",user_data->GetParser()->thistFile().c_str(),my_task_id);
      FILE *thist_fptr = fopen(my_thist_file_name,"w");
      fprintf(thist_fptr,"%s",results_str.c_str());
      fclose(thist_fptr);
    }

    MPI_Send(&worker_solution[0],
             msg_send_size,
             MPI_DOUBLE,
             MASTER_RANK,
             my_task_id,
             MPI_COMM_WORLD);


  } // end recv/send while loop
  //if(worker_solution != NULL) {
  //  delete [] worker_solution;
  //}

  if(afactor_mult != NULL) {
    delete [] afactor_mult;
  }
  if(user_data != NULL) {
    delete user_data;
  }
  // clean up cvode
  if(num_runs > 0) {
    delete [] worker_solution;
    N_VDestroy_Serial(cvode_state);
    CVodeFree(&cvode_memory);
    DestroyAuxCVodeStructs(aux_cvode);
  }
  return 0;
}

void PrintEstimateEndTime(const int num_done,
                          const int num_total,
                          const double start_time)
{
  double frac = (double)num_done/(double)num_total;
  double elapse = getHighResolutionTime()-start_time;

  printf("# [%5.1f%% complete] estimated remaining time is %10.2e seconds ...\n",
         frac*100.,elapse*(1.0-frac)/frac);
  fflush(stdout);
}

void SendWorkersKillTag(void)
{
  int num_threads;
  MPI_Comm_size(MPI_COMM_WORLD,&num_threads);

  for(int j=0; j<num_threads; ++j) {
    if(j != MASTER_RANK) {
      // syncronous send block the application until the receiver starts
      // to receive the message
      MPI_Ssend(0, 0, MPI_INT, j, KILL_TAG, MPI_COMM_WORLD);
    }
  }
  double start_time   = getHighResolutionTime();
  double current_time = start_time;
  while(current_time - start_time < SHUT_DOWN_TIME) {
    current_time = getHighResolutionTime();
  }
}

// --------------------------------------------------------------------------
//class FlushableIndexList
//{
// public:
//  FlushableIndexList(const size_t list_size);
//  ~FlushableIndexList() {};
//  size_t Flush(std::vector<size_t> *index_list);
//  // Flush returns the size of the index list that was flushed
//  bool Ready(const size_t index);
//  // returns true if the index was set to the ready for flush state, it will
//  // return false if the index was out-of-bounds
//
// private:
//  size_t list_size_;
//  size_t next_index_;
//  std::vector<bool> ready_;
//  std::vector<bool> flush_;
//};

FlushableIndexList::FlushableIndexList(const size_t list_size)
{
  list_size_ = list_size;
  next_index_ = 0;
  ready_.assign(list_size, false);
  flush_.assign(list_size, false);
}

size_t FlushableIndexList::Flush(std::vector<size_t> *index_list)
{
  size_t last_flush=0;
  bool any_flushed=false;
  index_list->clear();

  for(size_t j=next_index_; j<list_size_; ++j) {

    if(ready_[j]) {

      // set index as flushed and add to index list
      flush_[j] = true;
      index_list->push_back(j);
      last_flush = j;
      any_flushed = true;

    } else {

      break;
    }
  }
  if(any_flushed) {
    // if any elements were flushed, then set the next_index_ to check
    next_index_ = last_flush+1;
  }
  return index_list->size();
}

bool FlushableIndexList::Ready(const size_t index)
{
  if(index < list_size_) {
    ready_[index] = true;
    return true;
  }
  return false;
}

// class FlushableIndexList
//---------------------------------------------------------------------------
