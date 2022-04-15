#include <stdlib.h>
#include <stdio.h>

#include "AFactorIFP.h"
#include "idtControlParams.h"
#include "idtSolvers.h"
#include "utility_funcs.h"

#include "perturbAFactor_common.h"

using zerork::getHighResolutionTime;

int main(int argc, char *argv[])
{
  FILE *outFilePtr;
  int nIdtTemp,nSoln, nTask, remTask;
  double startTime, elapsedTime, finishTime;
  int *nthLargestId;
  double *idtCalcs;
  double *cpuTime;
  double *relSens;
  string headerInfo,sensInfo;
  int hrr1_id, hrr2_id;
  double delta_hrr1, delta_hrr2;

  checkCommandLine(argc,argv);

  idtControlParams idtCtrl(argv[1],1);

  nIdtTemp = idtCtrl.num_idt_temperatures_;
  nTask = idtCtrl.nRxn;
  nSoln = idtCtrl.num_results_;
  hrr1_id = nSoln-3;
  hrr2_id = nSoln-1;
  if(idtCtrl.doBothDir) {
    nSoln*=2;
    nTask*=2;
  }
  nTask++; // add an additional task for the unperturbed case
  remTask=nTask;
  nthLargestId = new int[idtCtrl.nRxn];
  idtCalcs = new double[nSoln*(idtCtrl.nRxn+1)];
  cpuTime  = new double[idtCtrl.nRxn+1];
  relSens  = new double[idtCtrl.nRxn*nIdtTemp];

  // write header to stdout
  getHeaderInfo(argc, argv, &idtCtrl, true, headerInfo);
  printf("%s",headerInfo.c_str()); fflush(stdout); 

  // open output file and write header to it
  outFilePtr = fopen(idtCtrl.outFile.c_str(),"w");
  if(outFilePtr == NULL) {
    printf("ERROR: can not open file %s for write operation\n",
           idtCtrl.outFile.c_str());
    exit(-1);
  }
  fprintf(outFilePtr,"%s",headerInfo.c_str()); fflush(outFilePtr);

  // print out the column headers
  getColHeaderInfo(&idtCtrl, headerInfo);
  printf("%s",headerInfo.c_str()); fflush(stdout); 

  startTime = getHighResolutionTime();
  solveIdtOriginal(&idtCtrl,
                   idtCalcs,
                   &cpuTime[0]);
  remTask--;
  elapsedTime = getHighResolutionTime()-startTime;
  finishTime = (double)nTask*cpuTime[0];

  if(idtCtrl.doBothDir) {
    for(int j=0; j<nSoln/2; j++) {
      // copy IDT times for unperturbed case to facilitate sensitivity
      // coefficient calculations
      idtCalcs[j+nSoln/2] = idtCalcs[j];
    }  
  }
  printf("%8d %12.3e %12.3e",0,cpuTime[0],finishTime);
  for(int j=0; j<nSoln; j++) {
    printf(" %20.11e",idtCalcs[j]);
  }
  printf("\n"); fflush(stdout);

  for(int j=1; j<=idtCtrl.nRxn; j++) {
    solveIdtPerturbRxn(j-1,
                       &idtCtrl,
                       &idtCalcs[j*nSoln],
                       &cpuTime[j]);
    remTask--;
    if(idtCtrl.doBothDir) {
      remTask--;
    }
    elapsedTime = getHighResolutionTime()-startTime;
    finishTime = (double)remTask*elapsedTime/(double)(nTask-remTask);
    printf("%8d %12.3e %12.3e",j,cpuTime[j],finishTime);
    for(int k=0; k<nSoln; k++) {
      printf(" %20.11e",idtCalcs[j*nSoln+k]);
    }
    printf("    %s\n",idtCtrl.mech->getReactionName(j-1)); fflush(stdout);
  } // for loop over all reactions
  
  //------------------------------------------------------------------------
  printf("# Writing sensitivity data to file %s\n",idtCtrl.outFile.c_str());

  calcRxnSensitivity(&idtCtrl,
                     &idtCalcs[0],
                     &idtCalcs[nSoln],
                     relSens,
                     nthLargestId);

  getColHeaderInfo_sensitivity(&idtCtrl,sensInfo);
  fprintf(outFilePtr,"%s",sensInfo.c_str()); fflush(outFilePtr);

  for(int j=0; j<idtCtrl.nRxn; j++) {

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
  fprintf(outFilePtr,"%s",headerInfo.c_str()); fflush(outFilePtr);
  for(int j=0; j<=idtCtrl.nRxn; j++) {

    fprintf(outFilePtr,"%8d %12.3e           --",j,cpuTime[j]);
    for(int k=0; k<nSoln; k++) {
      fprintf(outFilePtr," %20.11e",idtCalcs[j*nSoln+k]);
    }
    // report the delta of the two heat release rate calculations
    if(idtCtrl.doBothDir) {
      delta_hrr1 = idtCalcs[j*nSoln+hrr1_id] - 
        idtCalcs[j*nSoln+idtCtrl.num_results_+hrr1_id];
      delta_hrr2 = idtCalcs[j*nSoln+hrr2_id] -
        idtCalcs[j*nSoln+idtCtrl.num_results_+hrr2_id];
    } else {
      delta_hrr1 = idtCalcs[j*nSoln+hrr1_id]-idtCalcs[hrr1_id];
      delta_hrr2 = idtCalcs[j*nSoln+hrr2_id]-idtCalcs[hrr2_id];
    }
    fprintf(outFilePtr," %20.11e %20.11e", delta_hrr1, delta_hrr2);

    if(j==0) {
      fprintf(outFilePtr,"    orig (unperturbed)\n");
    }
    else {
      fprintf(outFilePtr,"    %s\n",idtCtrl.mech->getReactionName(j-1));
    }
    fflush(outFilePtr);
  }
  
  fclose(outFilePtr);
  delete [] nthLargestId;
  delete [] relSens;
  delete [] cpuTime;
  delete [] idtCalcs;
  return 0;
}
