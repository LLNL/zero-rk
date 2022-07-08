#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "AFactorIFP.h"
#include "idtControlParams.h"
#include "idtSolvers.h"
#include "utility_funcs.h"

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
  //printf("# DEBUG: before MPI_Init\n"); fflush(stdout);
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  //printf("# DEBUG: after MPI_Init, rank = %d\n", myRank); fflush(stdout);
   
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

  int hrr1_id, hrr2_id;
  double delta_hrr1, delta_hrr2;

  double startTime, elapsedTime,sumIdtTime;

  string headerInfo, sensInfo;
  
  FILE *outFilePtr;

  double *idtOrig;
  double cpuTimeOrig;

  MPI_Status status;
  double *msgResult;

  int *nthLargestId;
  double *idtWorker;
  double *cpuTimeWorker;
  double *relSens;
  
  startTime = getHighResolutionTime();
  idtControlParams idtCtrl(inpArgv[1],1); // 1 = print mech parser file

  // get worker pool
  MPI_Comm_size(MPI_COMM_WORLD,&nThread);
  nWorker=nThread-1;

  // initialization for the master thread
  // set key array data lengths
  nTask = idtCtrl.nRxn;
  nIdtTemp  = idtCtrl.num_idt_temperatures_;
  nSoln     = idtCtrl.num_results_;
  hrr1_id = nSoln-3;
  hrr2_id = nSoln-1;
  if(idtCtrl.doBothDir) {
    nSoln*=2;
  }
  msgRecvSize = nSoln+1; // add an extra datum for the cpuTime

  //printf("nIdtTemp = %d\n", nIdtTemp);
  //printf("nTask    = %d\n", nTask);
  //printf("nSoln    = %d\n", nSoln);
  // allocate data arrays
  idtOrig       = new double[nSoln];
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
  for(int j=1; j<=nWorker; j++) {
    MPI_Send(&taskId,         // message buffer (rxn idx to perturb)
             1,               // message buffer count
             MPI_INT,         // message buffer type
             j,               // process rank
             WORK_TAG,        // user-specified message tag
             MPI_COMM_WORLD); // default communicator for all threads
    taskId++;
  }
  
  // -------------------------------------------------------------------------
  // solve unperturbed mechanism
  //printf("# DEBUG: master, call solveIdtOriginal()\n"); fflush(stdout);
  solveIdtOriginal(&idtCtrl,
                   &idtOrig[0],
                   &cpuTimeOrig);
  //printf("# DEBUG: master, returned solveIdtOriginal()\n"); fflush(stdout);

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
    //printf("# DEBUG: Master received task id = %d\n",taskIdTag); 
    fflush(stdout);
    for(int j=0; j<nSoln; j++) {
      idtWorker[taskIdTag*nSoln+j]=msgResult[j];
    }
    cpuTimeWorker[taskIdTag]=msgResult[msgRecvSize-1];
    sumIdtTime+=cpuTimeWorker[taskIdTag];
     
    // assign the worker a new problem
    MPI_Send(&taskId,           // message buffer (rxn idx to perturb)
             1,                 // message buffer count
             MPI_INT,           // message buffer type
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
  calcRxnSensitivity(&idtCtrl,
                     &idtOrig[0],
                     &idtWorker[0],
                     relSens,
                     nthLargestId);
  elapsedTime=getHighResolutionTime()-startTime;
  fprintf(outFilePtr,"# Total elapsed time            [s]: %13.5e\n",
          elapsedTime);
  fprintf(outFilePtr,"# Sum of all IDT solution times [s]: %13.5e\n",
          sumIdtTime);
  fprintf(outFilePtr,"#----------------------------------------------------------------------------\n");

  printf("# Writing sensitivity data to file %s\n",idtCtrl.outFile.c_str());

  getColHeaderInfo_sensitivity(&idtCtrl,sensInfo);
  fprintf(outFilePtr,"%s",sensInfo.c_str()); fflush(outFilePtr);

  for(int j=0; j<nTask; j++) {

    fprintf(outFilePtr,"%8d",j+1);

    for(int k=0; k<nIdtTemp; k++) {
      fprintf(outFilePtr," %16.7e",relSens[j*nIdtTemp+k]);
    }
    fprintf(outFilePtr,"    %s\n",idtCtrl.mech->getReactionName(j));
    fflush(outFilePtr);
  }

  // add sorted data
  fprintf(outFilePtr,"\n\n");
  fprintf(outFilePtr,"# Relative sensitivities Srel are sorted by the maximum magnitude taken over\n");
  fprintf(outFilePtr,"# all %d IDT temperatures for each reaction.\n",nIdtTemp);
  fprintf(outFilePtr,"#------------------------------------------------------------------------------\n");
  fprintf(outFilePtr,"%s",sensInfo.c_str()); fflush(outFilePtr);
  for(int j=0; j<idtCtrl.nRxn; j++) {

    fprintf(outFilePtr,"%8d",nthLargestId[j]+1);

    for(int k=0; k<nIdtTemp; k++) {
      fprintf(outFilePtr," %16.7e",relSens[nthLargestId[j]*nIdtTemp+k]);
    }
    fprintf(outFilePtr,"    %s\n",
            idtCtrl.mech->getReactionName(nthLargestId[j]));
    fflush(outFilePtr);
  }
  
  // add raw data to file
  fprintf(outFilePtr,"\n\n");
  fprintf(outFilePtr,"# Raw ignition delay time data\n");
  getColHeaderInfo(&idtCtrl, headerInfo);
  fprintf(outFilePtr,"%s",headerInfo.c_str()); fflush(outFilePtr);
  fprintf(outFilePtr,"%8d %12.3e           --",0,cpuTimeOrig);
  for(int k=0; k<nSoln; k++) {
    fprintf(outFilePtr," %20.11e",idtOrig[k%idtCtrl.num_results_]);
  }
  fprintf(outFilePtr," %20.11e %20.11e    orig (unperturbed)\n", 0.0, 0.0);
        
  for(int j=0; j<nTask; j++) {

    fprintf(outFilePtr,"%8d %12.3e           --",j+1,cpuTimeWorker[j]);
    for(int k=0; k<nSoln; k++) {
      fprintf(outFilePtr," %20.11e",idtWorker[j*nSoln+k]);
    }
    // report the delta of the two heat release rate calculations
    if(idtCtrl.doBothDir) {
      delta_hrr1 = idtWorker[j*nSoln+hrr1_id] - 
        idtWorker[j*nSoln+idtCtrl.num_results_+hrr1_id];
      delta_hrr2 = idtWorker[j*nSoln+hrr2_id] -
        idtWorker[j*nSoln+idtCtrl.num_results_+hrr2_id];
    } else {
      delta_hrr1 = idtWorker[j*nSoln+hrr1_id]-idtOrig[hrr1_id];
      delta_hrr2 = idtWorker[j*nSoln+hrr2_id]-idtOrig[hrr2_id];
    }
    fprintf(outFilePtr," %20.11e %20.11e", delta_hrr1, delta_hrr2);

    fprintf(outFilePtr,"    %s\n",idtCtrl.mech->getReactionName(j));
    fflush(outFilePtr);
  }

  // free memory
  fclose(outFilePtr);
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
  
  MPI_Status status;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  checkCommandLine(inpArgc, inpArgv);
  //printf("Initialize idtControlParams on worker rank = %d\n",myRank);
  //fflush(stdout);
  idtControlParams idtCtrl(inpArgv[1],0); // 0 = don't print parser file
  //printf("Finished idtControlParams on worker rank = %d\n",myRank);  

  msgSendSize = idtCtrl.num_results_;
  if(idtCtrl.doBothDir) {
    msgSendSize*=2;
  }
  msgSendSize++; // add an extra datum for the cpuTime

  workerSoln = new double[msgSendSize];

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
    //printf("# DEBUG: worker = %d, task id = %d, call solveIdtPerturbRxn()\n",
    //       myRank, myTaskId);

    // perturb the reaction given by myTaskId
    solveIdtPerturbRxn(myTaskId,
                       &idtCtrl,
                       &workerSoln[0],
                       &workerSoln[msgSendSize-1]);
    //printf("# DEBUG: worker = %d, task id = %d, returned solveIdtPerturbRxn()\n",
    //       myRank, myTaskId);

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
    //printf("worker idt[0]: %.18g\n",workerSoln[0]);
  }

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
