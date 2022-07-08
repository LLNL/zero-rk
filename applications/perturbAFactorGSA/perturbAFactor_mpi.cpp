#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "GSA_AFactorIFP.h"
#include "idtControlParams.h"
#include "idtSolvers.h"
#include "utility_funcs.h"
#include "utilities/sequential_file_matrix.h"
#include "utilities.h"
#include "gsa_stats.h"

#include "perturbAFactor_common.h"

using zerork::getHighResolutionTime;

const int MASTER_RANK=0;
const int WORK_TAG=1;
const int KILL_TAG=2;
const int PRINT_FREQ=50;


void masterTasks(int inpArgc, char *inpArgv[]);
void workerTasks(int inpArgc, char *inpArgv[]);
void estimateEndTime(const int ndone,
                     const int ntotal,
                     const double startTime);

int main(int argc, char *argv[])
{
  int myRank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  
  if(myRank == MASTER_RANK) {
    masterTasks(argc, argv);
  }
  else {
    workerTasks(argc, argv);
  }
  MPI_Finalize(); 
  return 0;
}

void masterTasks(int inpArgc, char *inpArgv[])
{
  int nThread;
  int nWorker;
  int taskId;
  int taskIdTag;

  int nIdtTemp;
  int nSoln;
  int msgRecvSize;
  int nTask;

  double startTime, elapsedTime,sumIdtTime;

  string headerInfo, idtTaskInfo;
  
  FILE *outFilePtr, *checkFilePtr;

  double *idtOrig;
  double cpuTimeOrig;

  MPI_Status status;
  double *msgResult;

  int *nthLargestId;
  double *idtWorker;
  double *cpuTimeWorker;
  double *relSens;

  double *afactor_mult;
  
  startTime = getHighResolutionTime();
  idtControlParams idtCtrl(inpArgv[1],1); // 1 = print mech parser file
  zerork::utilities::SequentialFileMatrixRead gsa_matrix(idtCtrl.gsaMatrixFile.c_str());
  GSAStats stats(idtCtrl.nRxn);
  // check basic Sequential FileMatrix information
  // make sure the matrix size matches the mechanism
  if(gsa_matrix.num_columns() != idtCtrl.nRxn) {
    printf("ERROR: GSA Matrix file is incompatible with the mechanism.\n");
    printf("       The number of matrix columns must match the number of reactions.\n");
    printf("         GSA Matrix file %s number of columns   = %d\n",
           gsa_matrix.filename(),gsa_matrix.num_columns());
    printf("         Mechanism file %s  number of reactions = %d\n",
           idtCtrl.mechFile.c_str(),idtCtrl.nRxn);
    exit(-1);
  }
  // allocate array for the afactor_multiplier
  afactor_mult = new double[gsa_matrix.num_columns()];

  //// echo matrix as a check
  //printf("%5d %5d\n",gsa_matrix.num_rows(),gsa_matrix.num_columns());
  //for(int j=0; j<gsa_matrix.num_rows(); ++j) {
  //  gsa_matrix.GetNextRow(&afactor_mult[0]);
  //  for(int k=0; k<gsa_matrix.num_columns(); ++k) {
  //    printf("%5d %5d  %24.16e\n",j,k,afactor_mult[k]);
  //  } 
  //}
  
  // printf("#   Number of buffer rows = %d\n", gsa_matrix.max_buffer_rows());
  // printf("# GSA Matrix File Read = %s\n",gsa_matrix.filename());
  // printf("#   Matrix size (%d, %d)\n",
  //        gsa_matrix.num_rows(),
  //        gsa_matrix.num_columns());


  // get worker pool
  MPI_Comm_size(MPI_COMM_WORLD,&nThread);
  nWorker=nThread-1;
  if(nWorker < 1) {
    printf("ERROR: current version does not support MPI with %d worker threads.\n", nWorker);
    printf("       Exiting now.\n");
    exit(-1);
  }


  // initialization for the master thread
  // set key array data lengths
  nTask = gsa_matrix.num_rows();
  nSoln = nIdtTemp = idtCtrl.cvodeCtrl.nRoots;
  msgRecvSize = nSoln+1; // add an extra datum for the cpuTime

  //printf("nIdtTemp = %d\n", nIdtTemp);
  //printf("nTask    = %d\n", nTask);
  //printf("nSoln    = %d\n", nSoln);
  // allocate data arrays
  idtOrig       = new double[nIdtTemp];
  nthLargestId  = new int[nTask];
  idtWorker     = new double[nTask*nSoln];
  cpuTimeWorker = new double[nTask];
  relSens       = new double[nTask*nIdtTemp];
  msgResult     = new double[msgRecvSize];

  // -------------------------------------------------------------------------
  // assign the first batch of IDT problems to the workers
  if(nWorker > nTask) {
    printf("WARNING: the number of worker threads (%d) exceeds\n",
           nWorker);
    printf("         the number of AFactor perturbation tasks (%d).\n",
           nTask);
    nWorker=nTask;
  }
  taskId = 0;
  // CAUTION: assumes master rank is zero
  for(int j=1; j<=nWorker; j++) {

    gsa_matrix.GetNextRow(&afactor_mult[0]);
    stats.AddNextSample(&afactor_mult[0]);

    MPI_Send(&taskId,          // message buffer (task index)
             1,                // message buffer count
             MPI_INT,          // message buffer type
             j,                // process rank
             WORK_TAG,         // user-specified message tag
             MPI_COMM_WORLD);  // default communicator for all threads

    MPI_Send(&afactor_mult[0], // message buffer 
             idtCtrl.nRxn,     // message buffer count
             MPI_DOUBLE,       // message buffer type
             j,                // process rank
             WORK_TAG,         // user-specified message tag
             MPI_COMM_WORLD);  // default communicator for all threads


    taskId++;
  }
  
  // -------------------------------------------------------------------------
  // solve unperturbed mechanism
  solveIdtOriginal(&idtCtrl,
                   &idtOrig[0],
                   &cpuTimeOrig);

  sumIdtTime=cpuTimeOrig;

  // open output file and write header to it
  outFilePtr = fopen(idtCtrl.outFile.c_str(),"w");
  if(outFilePtr == NULL) {
    printf("ERROR: can not open file %s for write operation\n",
           idtCtrl.outFile.c_str());
    exit(-1);
  }
  getHeaderInfo(inpArgc, inpArgv, &idtCtrl, true, headerInfo);
  fprintf(outFilePtr,"%s",headerInfo.c_str()); fflush(outFilePtr);

  // -------------------------------------------------------------------------
  // receive the results and send out the remaining tasks
  while(taskId < nTask) {
    // recieve result from any worker
    MPI_Recv(&msgResult[0],     // results buffer
             msgRecvSize,       // size of buffer
             MPI_DOUBLE,        // buffer type
             MPI_ANY_SOURCE,    // receive msgResults from any thread
             MPI_ANY_TAG,       // receive msgResults with any tags
             MPI_COMM_WORLD,    // default communicator for all threads
             &status);          // info about the communication

    // record the results for the task id
    taskIdTag = status.MPI_TAG; // workers send the taskId in the MPI_TAG

    for(int j=0; j<nSoln; j++) {
      idtWorker[taskIdTag*nSoln+j]=msgResult[j];
    }
    cpuTimeWorker[taskIdTag]=msgResult[msgRecvSize-1];
    sumIdtTime+=cpuTimeWorker[taskIdTag];

    // get the next set of A-Factor pertubation multipliers
    gsa_matrix.GetNextRow(&afactor_mult[0]);
    stats.AddNextSample(&afactor_mult[0]);
     
    // assign the worker a new problem
    MPI_Send(&taskId,           // message buffer (rxn idx to perturb)
             1,                 // message buffer count
             MPI_INT,           // message buffer type
             status.MPI_SOURCE, // process rank that just finished
             WORK_TAG,          // user-specified message tag
             MPI_COMM_WORLD);   // default communicator for all threads
 

    MPI_Send(&afactor_mult[0],  // message buffer 
             idtCtrl.nRxn,      // message buffer count
             MPI_DOUBLE,        // message buffer type
             status.MPI_SOURCE, // process rank that just finished
             WORK_TAG,          // user-specified message tag
             MPI_COMM_WORLD);   // default communicator for all threads

    taskId++;
    if(taskId%PRINT_FREQ == 0) {
      estimateEndTime(taskId,nTask,startTime);
    }
  }      
  // -------------------------------------------------------------------------
  // All tasks have been assigned, and all workers (rank 1,...,nWorker) 
  // received at least one task.  Now collect the last results from each
  // worker.
  for(int j=1; j<=nWorker; j++) {
    MPI_Recv(&msgResult[0],     // results buffer
             msgRecvSize,       // size of buffer
             MPI_DOUBLE,        // buffer type
             MPI_ANY_SOURCE,    // receive msgResults from any thread
             MPI_ANY_TAG,       // receive msgResults with any tags
             MPI_COMM_WORLD,    // default communicator for all threads
             &status);          // info about the communication

    // record the results for the task id
    taskIdTag = status.MPI_TAG; // workers send the taskId in the MPI_TAG

    for(int j=0; j<nSoln; j++) {
      idtWorker[taskIdTag*nSoln+j]=msgResult[j];
    }
    cpuTimeWorker[taskIdTag]=msgResult[msgRecvSize-1];
    sumIdtTime+=cpuTimeWorker[taskIdTag];
  }    

  // -------------------------------------------------------------------------
  // no more work - send out the worker KILL_TAG with an empty message
  for(int j=1; j<=nWorker; j++) {
    MPI_Send(0, 0, MPI_INT, j, KILL_TAG, MPI_COMM_WORLD);
  }

  // -------------------------------------------------------------------------
  // post-process the results collected

  elapsedTime=getHighResolutionTime()-startTime;
  fprintf(outFilePtr,"# Total elapsed time            [s]: %13.5e\n",
          elapsedTime);
  fprintf(outFilePtr,"# Sum of all IDT solution times [s]: %13.5e\n",
          sumIdtTime);
  fprintf(outFilePtr,"#----------------------------------------------------------------------------\n");

  printf("# Writing ignition delay data to file %s\n",idtCtrl.outFile.c_str());

  // write the data file header
  getColHeaderInfo_IdtTask(&idtCtrl,idtTaskInfo);
  fprintf(outFilePtr,"%s",idtTaskInfo.c_str()); fflush(outFilePtr);

  // report the unperturbed IDT calculation as task 0
  fprintf(outFilePtr,"%8d %16.7e",0,cpuTimeOrig);
  for(int k=0; k<nSoln; ++k) {
    fprintf(outFilePtr," %16.7e",idtOrig[k]);
  }
  fprintf(outFilePtr,"\n");

  for(int j=0; j<nTask; j++) {

    fprintf(outFilePtr,"%8d %16.7e",j+1,cpuTimeWorker[j]);

    for(int k=0; k<nSoln; k++) {
      fprintf(outFilePtr," %16.7e",idtWorker[j*nSoln+k]);
    }
    fprintf(outFilePtr,"\n");
    fflush(outFilePtr);
  }

  // reaction check for bounds min, max A-Factor and estimate
  // of normal parameters for log(afactor_mult)
  checkFilePtr = fopen(idtCtrl.checkFile.c_str(),"w");
  if(checkFilePtr == NULL) {
    printf("ERROR: can not open the stat check file %s for write operation\n",
           idtCtrl.checkFile.c_str());
    exit(-1);
  }
  fprintf(checkFilePtr,
          "# Basic distribution statistics for each reaction's perturbation multiplier.\n");
  fprintf(checkFilePtr,
          "# Here 'a' is the perturbation multiplier, so a=1 is unperturbed.\n");
  fprintf(checkFilePtr,"# Note reaction id is index-0.\n");
  fprintf(checkFilePtr,"# Number of samples   = %d (GSA matrix rows)\n",
          stats.num_samples());
  fprintf(checkFilePtr,"# Number of reactions = %d (GSA matrix columns)\n",
          stats.num_dimensions());
  fprintf(checkFilePtr,"# -----------------------------------------------------------------------------\n");
  fprintf(checkFilePtr,"#   id  reaction name                            min(a)         max(a)   mean[log(a)]  stdev[log(a)]\n");
  for(int j=0; j<idtCtrl.nRxn; ++j) {
    // note that '-' indicates left-justify for value
    fprintf(checkFilePtr,"%6d  %-32s  %13.4e  %13.4e  %13.4e  %13.4e\n",
            j,
            idtCtrl.mech->getReactionName(j),
            stats.MinValue(j),
            stats.MaxValue(j),
            stats.LogMean(j),
            sqrt(stats.LogVariance(j)));
  }



  // free memory
  fclose(outFilePtr);
  fclose(checkFilePtr);

  delete [] afactor_mult;  
  delete [] msgResult;
  delete [] relSens;
  delete [] cpuTimeWorker;
  delete [] idtWorker;
  delete [] nthLargestId;
  delete [] idtOrig;  

}
void workerTasks(int inpArgc, char *inpArgv[])
{
  int myRank;
  int myTaskId;
  int idtRuns=0;
  int msgSendSize;
  double *workerSoln;

  double *afactor_mult;
  
  MPI_Status status;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  checkCommandLine(inpArgc, inpArgv);
  //printf("Initialize idtControlParams on worker rank = %d\n",myRank);
  //fflush(stdout);
  idtControlParams idtCtrl(inpArgv[1],0); // 0 = don't print parser file
  //printf("Finished idtControlParams on worker rank = %d\n",myRank);  

  msgSendSize = idtCtrl.cvodeCtrl.nRoots;
  msgSendSize++; // add an extra datum for the cpuTime

  workerSoln = new double[msgSendSize];
  afactor_mult = new double[idtCtrl.nRxn];

  while(1) {
    MPI_Recv(&myTaskId,       // buffer - rxn id to perturb
             1,               // buffer size
             MPI_INT,         // buffer type
             MASTER_RANK,     // only receive messages from the MASTER
             MPI_ANY_TAG,     // receive messages with any tag
             MPI_COMM_WORLD,  // use the default communicator with all threads
             &status);        // info about the message

    if(status.MPI_TAG == KILL_TAG) {
      break; // leave the loop when the master thread send the KILL_TAG
    }
    MPI_Recv(&afactor_mult[0], // message buffer 
             idtCtrl.nRxn,     // message buffer count
             MPI_DOUBLE,       // message buffer type
             MASTER_RANK,      // only receive messages from the MASTER
             WORK_TAG,         // user-specified message tag
             MPI_COMM_WORLD,
             &status);  // default communicator for all threads

    // solve the ignition delay for the conditions
    solveIdtGsaPerturbRxn(afactor_mult,
                          &idtCtrl,
                          &workerSoln[0],
                          &workerSoln[msgSendSize-1]);
    idtRuns++;

    MPI_Send(&workerSoln[0],
             msgSendSize,
             MPI_DOUBLE,
             MASTER_RANK,
             myTaskId,
             MPI_COMM_WORLD);


  } // end recv/send while loop 

  if(idtRuns == 0) {
    // call the IDT solver at least once to avoid the segfault in the
    // destructor related to the SuperLU structrures.
    solveIdtOriginal(&idtCtrl,
                     &workerSoln[0],
                     &workerSoln[msgSendSize-1]);
    printf("INFO: worker rank %d performed no IDT calculations\n",myRank);
    printf("      computing one IDT calculations to avoid destructor\n");
    printf("      segmentation faults.\n");
    printf("      worker idt[0] =  %.18g for dT[0] = %.18g\n",
           workerSoln[0],
           idtCtrl.odeUserParams.redTempRoot[0]*idtCtrl.refTemp -
           idtCtrl.initTemp);
    fflush(stdout);
  }
  delete [] afactor_mult;
  delete [] workerSoln;
}

void estimateEndTime(const int ndone,
                     const int ntotal,
                     const double startTime)
{
  double frac = (double)ndone/(double)ntotal;
  double elapse = getHighResolutionTime()-startTime;
  
  printf("[%4.1f%% complete] estimated remaining time is %10.2e seconds ...\n",
         frac*100.,elapse*(1.0-frac)/frac);
  fflush(stdout);

}
