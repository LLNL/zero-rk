///////////////////////////////////////////////////////////////////////
//
//     This demonstration program builds an object representing a
//     reacting gas mixture, and uses it to compute thermodynamic
//     properties, chemical equilibrium, and transport properties.
//
///////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const int MAX_SPECNAME_LEN=256;
const int MAX_FLOAT_LEN=256;
const int MAX_THIST_LINE=20000;

//#include <culapack.h>

#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
//#include <cvode/cvode_spgmr.h>      // prototypes & constants for CVSPGMR
#include <cvode/cvode_dense.h>     // prototype for CVLapackDense
#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros

#include "CKconverter/CKReader.h"
#include "zerork/constants.h"

#include "ode_funcs.h"
#include "sweep_util.h"
#include "utility_funcs.h"


int setMoleFracFromFile(zerork::mechanism &mechInp, const char *fileName,
			double moleFrac[]);

static int check_flag(void *flagvalue, const char *funcname, int opt);

void getTimeHistLine_full(void *cv_mem,
                          const double currTime, 
                          const double *sysState,
                          const idt_sweep_params *idt_ctrl, 
                          const ODEParams_CMT *cvp,
                          const double wallTime,
                          char *thistLine);

void getTimeHistHeader_full(const idt_sweep_params *idt_ctrl,
                            const ODEParams_CMT *cvp,
                            char *thistHead);
int getIdtLine(const double initMoleFrac[],
               const idt_sweep_params *idt_ctrl,
               const double idt,
               const double runTime,
               const int nCvodeFails,
               const int maxBuffer,
               char *idtLine);
int getIdtHeader(const idt_sweep_params *idt_ctrl,
                 const ODEParams_CMT *cvp,
                 const int maxBuffer,
                 char *idtHeader);
                

int getSimHeader_full(const int inp_argc,
                      char **inp_argv,
                      const idt_sweep_params *idt_ctrl,
                      const ODEParams_CMT *cvp,
                      const int max,
                      char *simHead);


void cvReactor(int inp_argc, char **inp_argv)
{
  double idt,tmax,dtprint;
  //int isIdtFound;
  double *avgTime, *minTime, *maxTime;
  double nRunsPerThresh;
  double concentration_ref;
  double temperature_roots[1];

  FILE *idtFilePtr,*thistFilePtr;
  int j;
  int nSpc,nState;
  double *moleFracInit;
  double *moleFracCurr;
  double *concentration_p;
  double tnext,tcurr;
  // cvode variables
  N_Vector systemState;
  N_Vector system_derivative;
  ODEParams_CMT *systemParams;
  void *cvode_mem;
  int flag;

  char thistLine[MAX_THIST_LINE];
  int did_cvodeFail;
  int isBadStep;

  // timing data
  double startTime,stopTime, simTime;
//  double setupZeroRK,setupCVode,setupJterm;

  // sweep controls
  idt_sweep_params *idt_ctrl;

  if(inp_argc != 2)
    {
      printf("ERROR: incorrect command line usage\n");
      printf("       use instead %s <idt sweep input>\n",inp_argv[0]);
      fflush(stdout); exit(-1);
    }

  startTime=getHighResolutionTime();
  idt_ctrl=new idt_sweep_params(inp_argv[1]);
  stopTime=getHighResolutionTime();
//  setupZeroRK=stopTime-startTime;

  nSpc=idt_ctrl->getNumSpecies();
  nState=idt_ctrl->getNumSpecies()+2;
  //nStep=idt_ctrl->getNumSteps();

  avgTime=(double *)malloc(sizeof(double)*idt_ctrl->getNumThreshRuns());
  minTime=(double *)malloc(sizeof(double)*idt_ctrl->getNumThreshRuns());
  maxTime=(double *)malloc(sizeof(double)*idt_ctrl->getNumThreshRuns());

  for(j=0; j<idt_ctrl->getNumThreshRuns(); j++)
    {
      avgTime[j]=0.0;
      minTime[j]=1.0e300;
      maxTime[j]=-1.0e300;
    }

  systemState = N_VNew_Serial(nState);
  system_derivative = N_VNew_Serial(nState);
  concentration_p=NV_DATA_S(systemState); // caution: assumes realtype == double
  moleFracCurr = (double *)malloc(sizeof(double)*nSpc);
  moleFracInit = (double *)malloc(sizeof(double)*nSpc);

  
  // set up the system parameters
  concentration_ref = idt_ctrl->getInitPres()/
    (idt_ctrl->getMechPtr()->getGasConstant()*idt_ctrl->getInitTemp());

  temperature_roots[0] = idt_ctrl->getInitTemp()+idt_ctrl->getDeltaTign();
  systemParams = AllocateODEParams_CMT(idt_ctrl->getMechPtr(),
                                       concentration_ref,
                                       idt_ctrl->getRefTemp(),
                                       1,
                                       temperature_roots);
                                       
  // reset the system state
  idt_ctrl->getInitMoleFrac(concentration_p);
  NV_Ith_S(systemState,nSpc)   = 1.0; //
  NV_Ith_S(systemState,nSpc+1) = idt_ctrl->getInitTemp()/
    idt_ctrl->getRefTemp();
  ConstantVolumeWSR_CMT(0.0,systemState,system_derivative,systemParams);
  //exit(-1);
  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
  startTime=getHighResolutionTime();
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) exit(-1);
  //systemParam.cvodeMemPtr=cvode_mem;
 
  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag=CVodeInit(cvode_mem, ConstantVolumeWSR_CMT, 0.0, systemState);
  if (check_flag(&flag, "CVodeInit", 1)) exit(-1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  flag = CVodeSStolerances(cvode_mem,idt_ctrl->getRTol(),idt_ctrl->getATol());
  if (check_flag(&flag, "CVodeSStolerances", 1)) exit(-1);

  /* Call CVodeRootInit to specify the root function with 1 component */
  flag = CVodeRootInit(cvode_mem, 1, TempRootFunc_CMT);
  //flag = CVodeRootInit(cvode_mem, 2, tempRoot2);
  if (check_flag(&flag, "CVodeRootInit", 1)) exit(-1);

  /* Set the pointer to user-defined data */
  flag = CVodeSetUserData(cvode_mem, systemParams);
  if(check_flag(&flag, "CVodeSetUserData", 1)) exit(-1);

   /* Set the pointer to user-defined data */
  flag = CVDense(cvode_mem, nState);
  if(check_flag(&flag, "CVDense", 1)) exit(-1);
 

  // /* Call CVSpgmr to specify the linear solver CVSPGMR 
  //    with left preconditioning and the maximum Krylov dimension maxl */
  // flag = CVSpgmr(cvode_mem, PREC_LEFT, idt_ctrl->getKrylovDim());
  // if(check_flag(&flag, "CVSpgmr", 1)) exit(-1);

  // flag = CVSpilsSetGSType(cvode_mem, MODIFIED_GS);
  // if(check_flag(&flag, "CVSpilsSetGSType", 1)) exit(-1);

  // /* Set preconditioner setup and solve routines Precond and PSolve,
  //    and the pointer to the user-defined block data */
  // flag = CVSpilsSetPreconditioner(cvode_mem, jac_full_prec_setup,
  // 				  jac_full_prec_solveV3);
  // if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) exit(-1);

//  flag = CVSpilsSetJacTimesVecFn(cvode_mem, sparse_jac_v);
//  if(check_flag(&flag, "CVSpilsSetJacTimesVecFn", 1)) exit(-1);

  /* Set the maximum number of internal steps per CVode call and the maximum
   * allowable internal steps. */
  flag = CVodeSetMaxNumSteps(cvode_mem, idt_ctrl->getMaxInternalSteps());
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) exit(-1);
  flag = CVodeSetMaxStep(cvode_mem, idt_ctrl->getMaxInternalDt());
  if (check_flag(&flag, "CVodeSetMaxStep", 1)) exit(-1);

  flag = CVodeSetNonlinConvCoef(cvode_mem, idt_ctrl->getNlConvCoeff()); // Default [0.1]
  //flag = CVSpilsSetEpsLin(cvode_mem, idt_ctrl->getEpsLin());    // Default [0.05]

  stopTime=getHighResolutionTime();
  // setupCVode=stopTime-startTime;
  //nTrackSpec = idt_ctrl->getNumTrackSpecies();

  // ready for integration

  // open ignition delay and time history files for write
  idtFilePtr=fopen(idt_ctrl->getIdtFileName(),"w");
  if(idtFilePtr==NULL)
    {
      printf("ERROR: could not open output file %s for write\n",
	     idt_ctrl->getIdtFileName());
      exit(-1);
    }
  thistFilePtr=fopen(idt_ctrl->getTHistFileName(),"w");
  if(thistFilePtr==NULL)
    {
      printf("ERROR: could not open output file %s for write\n",
	     idt_ctrl->getTHistFileName());
      exit(-1);
    }
  // create the simulation header stored in a character array
  getSimHeader_full(inp_argc,
                    &inp_argv[0],
                    idt_ctrl,
                    systemParams,
                    MAX_THIST_LINE,
                    thistLine);
  // write the simulation  header to stdout, idt file and thist file
  printf("%s",thistLine);
  fflush(stdout);
  fprintf(thistFilePtr,"%s",thistLine); 
  fflush(thistFilePtr);
  fprintf(idtFilePtr,"%s",thistLine); 
  fflush(idtFilePtr);

  // create the time history header stored in a character array
  getTimeHistHeader_full(idt_ctrl, systemParams, thistLine); 
  // write the time history header to stdout and the time history file 
  printf("%s",thistLine);
  fflush(stdout);
  fprintf(thistFilePtr,"%s",thistLine); 
  fflush(thistFilePtr);

  // create the idt file header stored in a character array
  getIdtHeader(idt_ctrl, systemParams, MAX_THIST_LINE, thistLine);
  // write the idt header to the idt file
  fprintf(idtFilePtr,"%s",thistLine);
  fflush(idtFilePtr);

  tmax    = idt_ctrl->getStopTime();
  dtprint = idt_ctrl->getPrintTime();

  for(j=0; j<idt_ctrl->getRunTotal(); j++) {

    // reset the system parameters
    concentration_ref = idt_ctrl->getInitPres()/
      (idt_ctrl->getMechPtr()->getGasConstant()*idt_ctrl->getInitTemp());
    ResetODEParams_CMT(concentration_ref,
                       systemParams);
    // initialize the run counters
    //isIdtFound=0;
    did_cvodeFail=0;

    // reset the ignition delay threshold
    systemParams->temperature_roots[0] =
      idt_ctrl->getDeltaTign()+idt_ctrl->getInitTemp();

    // reset the system state
    idt_ctrl->getInitMoleFrac(concentration_p);
    NV_Ith_S(systemState,nSpc)   = 1.0; //
    NV_Ith_S(systemState,nSpc+1) = idt_ctrl->getInitTemp()/
      idt_ctrl->getRefTemp();

    // reset the time
    tcurr=0.0;
    tnext=dtprint;
    idt=tmax;
 
    // reinitialize cvode
    flag = CVodeReInit(cvode_mem, tcurr, systemState);
    isBadStep = 0;

    getTimeHistLine_full(cvode_mem,
                         tcurr,
                         NV_DATA_S(systemState),
                         idt_ctrl,
                         systemParams,
                         0.0,
                         thistLine);  
    printf("%s",thistLine);
    fflush(stdout);
    fprintf(thistFilePtr,"%s",thistLine);
    fflush(thistFilePtr);
    startTime=getHighResolutionTime();
 
    while(tcurr<tmax) {

      //printf("tcurr = %.18g\n",tcurr); fflush(stdout);    

      if(idt_ctrl->oneStep()) {
        flag = CVode(cvode_mem, tmax, systemState, &tcurr, CV_ONE_STEP);
      } else {
  	flag = CVode(cvode_mem, tnext, systemState, &tcurr, CV_NORMAL);
      }

      if (check_flag(&flag, "CVODE ERROR", 1)) {

        isBadStep = 1;
  	did_cvodeFail++;
        if(did_cvodeFail <= idt_ctrl->getMaxPrimaryCvodeFails()) {

  	  printf("WARNING: attempting CVodeReInit at t=%.18g [s]\n",tcurr);
  	  printf("         after cvode failure count = %d\n",
          did_cvodeFail);
	  
          flag = CVodeReInit(cvode_mem, tcurr, systemState);
        } else {

          printf("ERROR: idt calculation failed to recover after\n");
          printf("       %d re-initializations.  Halting Now!\n",
                 did_cvodeFail-1);
          exit(-1);
        }
      } // end if(check_flag(&flag, "CVODE ERROR", 1))		  

    getTimeHistLine_full(cvode_mem,
                         tcurr,
                         NV_DATA_S(systemState),
                         idt_ctrl,
                         systemParams,
                         getHighResolutionTime()-startTime,
                         thistLine);  

     printf("%s",thistLine);
     fflush(stdout);
     fprintf(thistFilePtr,"%s",thistLine);
     fflush(thistFilePtr);

      if(flag == CV_ROOT_RETURN) {
        idt=tcurr;
        //isIdtFound=1;
        if(idt_ctrl->continueAfterIDT() == 0) {
          tcurr=2.0*tmax; // advance time past the stopping criteria
         }
      }

  	  else if(!idt_ctrl->oneStep() && isBadStep==0)
  	    {tnext+=dtprint;}	  
  	}

    stopTime=getHighResolutionTime();
    simTime=stopTime-startTime;
    // add new-line to separate the time histories
    printf("\n");
    fprintf(thistFilePtr,"\n");

    idt_ctrl->getInitMoleFrac(moleFracInit);
    getIdtLine(moleFracInit,idt_ctrl,idt,simTime,did_cvodeFail,
               MAX_THIST_LINE,thistLine);
    fprintf(idtFilePtr,"%s",thistLine);
    fflush(idtFilePtr);

    avgTime[idt_ctrl->getThreshId()]+=simTime;
    if(simTime > maxTime[idt_ctrl->getThreshId()]) {
      maxTime[idt_ctrl->getThreshId()]=simTime;
    }
    if(simTime < minTime[idt_ctrl->getThreshId()]) {
      minTime[idt_ctrl->getThreshId()]=simTime;
    }

    idt_ctrl->incrementRunId();
  } // end for j-loop

  // report the run statistics
  fprintf(idtFilePtr,"\n\n");
  fprintf(idtFilePtr,"# -----------------------------------------------------------------------------\n");
  fprintf(idtFilePtr,"# Column 1: [#] number of species\n");
  fprintf(idtFilePtr,"# Column 2: [#] length of ODE state vector\n");
  fprintf(idtFilePtr,"# Column 3: [-] preconditioner threshold\n");
  fprintf(idtFilePtr,"# Column 4: [s] average simulation time\n");
  fprintf(idtFilePtr,"# Column 5: [s] minimum simulation time\n");
  fprintf(idtFilePtr,"# Column 6: [s] maximum simulation time\n");
  nRunsPerThresh=(double)(idt_ctrl->getRunTotal()/
			  idt_ctrl->getNumThreshRuns());

  for(j=0; j<idt_ctrl->getNumThreshRuns(); j++) {
    avgTime[j]/=nRunsPerThresh;
    fprintf(idtFilePtr,"%6d  %6d  %14.7e  %14.7e  %14.7e  %14.7e\n",
            nSpc,nState,idt_ctrl->getThreshFromId(j),
            avgTime[j],minTime[j],maxTime[j]);
  }
  fflush(idtFilePtr);

  // clean up procedures
  fclose(idtFilePtr);
  fclose(thistFilePtr);
 // Free ode system parameter memory
  // free(systemParam.molWt);
  // free(systemParam.invMolWt);
  // free(systemParam.netProd);
  // free(systemParam.fwdROP);
  // free(systemParam.createRate);
  // free(systemParam.destroyRate);
  // free(systemParam.conc);
  // free(systemParam.Energy);
  // free(systemParam.CvMass);

  // free_Jsparse(systemParam.sparseMtx);
  free(avgTime);
  free(minTime);
  free(maxTime);
  free(moleFracCurr);
  free(moleFracInit);
  N_VDestroy_Serial(systemState);
  N_VDestroy_Serial(system_derivative);

  FreeODEParams_CMT(systemParams);
  // Free integrator memory
  CVodeFree(&cvode_mem);

//  printf("Begin.......delete of idt_sweep_params class\n"); fflush(stdout);
  delete idt_ctrl;
//  printf("Finished....delete of idt_sweep_params class\n"); fflush(stdout);
}
     
int main(int argc, char **argv)
{
  cvReactor(argc, argv);
  exit(0);
}

// getTimeHistLine_full(...)
//
// The purpose of this function is to construct a data line stored in a
// character pointer *thistLine.  The data line is considered 'full' because
// it contains all the ODE solver statistics and complete timing information.
// The data recorded per line is as follows: run id, ode time, temperature,
// pressure, density,  mole fraction of tracked species (i.e. those specified
// in the input file), ODE solver procedure statistics, SuperLU statistics,
// solver timings.
//
// Note that changes to getTimeHistLine_full(...) should be matched with
// getTimeHistHeader_full(...) to ensure that the labels of the data columns
// are consistent.
//
// TODO: Currently, there is no check of the array size so it is possible to
// overwrite the data, which needs to be fixed.
void getTimeHistLine_full(void *cv_mem,
                          const double currTime, 
                          const double *sysState,
                          const idt_sweep_params *idt_ctrl, 
                          const ODEParams_CMT *cvp,
                          const double wall_time,
                          char *thistLine)
{
  int j;
  int nSpc   = cvp->mechp->getNumSpecies();
  int nTrack = idt_ctrl->getNumTrackSpecies();
  int strPos,strLen;
  double currPres,currTemp,currDens, initTemp;
  double *currMoleFrac;
  long int nsteps,netfails,nniters,nncfails;

  double molWt, cvMixPerMass, intEnergyMixPerMass, intEnergyMixPerMassInit;

  nsteps=netfails=nniters=nncfails=0;

  CVodeGetNumSteps(cv_mem, &nsteps);
  CVodeGetNumErrTestFails(cv_mem, &netfails);
  CVodeGetNumNonlinSolvIters(cv_mem, &nniters);
  CVodeGetNumNonlinSolvConvFails(cv_mem, &nncfails);

  currMoleFrac = new double[nSpc];

  initTemp = idt_ctrl->getInitTemp();
  currTemp = sysState[nSpc+1]*idt_ctrl->getRefTemp();
  currPres = sysState[nSpc]*         // p = C_{mix}*R_u*T
    cvp->concentration_ref*
    cvp->mechp->getGasConstant()*
    currTemp;
 
  molWt=cvp->mechp->getMolWtMixFromY(&sysState[0]);
  currDens = molWt*sysState[nSpc]*cvp->concentration_ref;
  // currMoleFrac[] used as temporary storage for species specific heat
  cvMixPerMass=cvp->mechp->getMassCvFromTY(currTemp,
                                          &sysState[0],
                                          &currMoleFrac[0]);
  // currMoleFrac[] used as temporary storage for species internal energy
  intEnergyMixPerMass=cvp->mechp->getMassIntEnergyFromTY(currTemp,
                                                        &sysState[0],
                                                        &currMoleFrac[0]);
  // currMoleFrac[] used as temporary storage for species internal energy
  intEnergyMixPerMassInit=cvp->mechp->getMassIntEnergyFromTY(initTemp,
                                                            &sysState[0],
                                                            &currMoleFrac[0]);
  // currMoleFrac[] is recalculated with mole fraction
  cvp->mechp->getXfromY(&sysState[0],currMoleFrac);

  // sprintf return the number of characters written excluding the '\0'
  // string termination character
  strLen=sprintf(thistLine,
		 "%8d  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
		 idt_ctrl->getRunId(),
		 currTime,
		 currTemp,
		 currPres,
		 currDens,
                 molWt,
                 cvMixPerMass,
                 intEnergyMixPerMass-intEnergyMixPerMassInit);
  strPos=strLen;
  for(j=0; j<nTrack; j++)
    {
      strLen=sprintf(&thistLine[strPos],
              "  %14.7e",
              currMoleFrac[idt_ctrl->getTrackSpeciesId(j)]);
      strPos+=strLen;
    }

  strLen=sprintf(&thistLine[strPos],
                 "  %6ld  %6ld  %6ld  %6ld",
                 nsteps,
                 netfails,
                 nniters,
                 nncfails);
  strPos+=strLen;
  strLen=sprintf(&thistLine[strPos],
                 "  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e\n",
                 cvp->timer.jacobian_time,
                 cvp->timer.lu_factor_time,
                 cvp->timer.lu_solve_time,
                 cvp->timer.function_time,
                 wall_time);

  delete [] currMoleFrac;
}

// getTimeHistHeader_full(...)
//
// The purpose of this function is to construct a descriptive header stored in 
// a character pointer *thistHead.  The header provides labels of each column
// written by getTimeHistLine_full(...).  The header and data line are
// considered 'full' because they contain all the ODE solver statistics and
// complete timing information. The data recorded per line is as follows:
// run id, ode time, temperature, pressure, density,  mole fraction of tracked
// species (i.e. those specified in the input file), ODE solver procedure
// statistics, SuperLU statistics, solver timings.
//
// Note that changes to getTimeHistLine_full(...) should be matched with
// getTimeHistHeader_full(...) to ensure that the labels of the data columns
// are consistent.
//
// TODO: Currently, there is no check of the array size so it is possible to
// overwrite the data, which needs to be fixed.
void getTimeHistHeader_full(const idt_sweep_params *idt_ctrl,
                            const ODEParams_CMT *cvp,
                            char *thistHead)
{
  int j;
  //int nSpc   = cvp->mechp->getNumSpecies();
  int nTrack = idt_ctrl->getNumTrackSpecies();
  int strPos,strLen, colNum;

  strPos=0;
  strLen=sprintf(&thistHead[strPos],
                 "%s%s%s%s%s%s%s%s",
                 "# Column  1: [#]       run id\n",
                 "# Column  2: [s]       simulated ODE time\n",
                 "# Column  3: [K]       temperature\n",
                 "# Column  4: [Pa]      pressure\n",
                 "# Column  5: [kg/m^3]  density\n",
                 "# Column  6: [kg/kmol] mixture molecular weight\n",
                 "# Column  7: [J/kg/K]  mixture specific heat (cosnt vol)\n",
                 "# Column  8: [J/kg]    energy released\n");
  strPos+=strLen;
  colNum=9;
  for(j=0; j<nTrack; j++) {
    strLen=sprintf(&thistHead[strPos],
                   "# Column %2d: [-]      mole fraction of %s\n",
                   colNum++, // increment after printf
                   cvp->mechp->getSpeciesName(idt_ctrl->getTrackSpeciesId(j)));
    strPos+=strLen;
  }
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [#]      CVodeGetNumSteps(...)\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [#]      CVodeGetNumErrTestFails(...)\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [#]      CVodeGetNumNonlinSolvIters(...)\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [#]      CVodeGetNumNonlinSolvConvFails(...)\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [s]      wall clock time - Jacobian setup\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [s]      wall clock time - LU factorization\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [s]      wall clock time - LU backsolve\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [s]      wall clock time - RHS evaluation\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [s]      wall clock time - total simulation\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# -----------------------------------------------------------------------------\n");
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# run id        sim time            temp           press            dens          mol wt              cv       E release");
  strPos+=strLen;
  for(j=0; j<nTrack; j++) {
    strLen=sprintf(&thistHead[strPos],
                   " mlfrc%10s",
                   cvp->mechp->getSpeciesName(idt_ctrl->getTrackSpeciesId(j)));
    strPos+=strLen;
  }
  strLen=sprintf(&thistHead[strPos],
                 "   nstep   nertf  nnlsit  nnlcvf    jac setup tm       LU fac tm     LU solve tm     RHS func tm     Total sim tm\n");
  strPos+=strLen;
}

//
int getSimHeader_full(const int inp_argc,
                      char **inp_argv,
                      const idt_sweep_params *idt_ctrl,
                      const ODEParams_CMT *cvp,
                      const int maxLen,
                      char *simHead)
{
  int j,spcId;
  int nSpc   = cvp->mechp->getNumSpecies();
  //int nSpcTermsJac = cvp->sparseMtx->nNonZero-2*nSpc+1;
  //int nTrack = idt_ctrl->getNumTrackSpecies();
  int strPos,strLen;

  // 1. construct the command line
  strPos=0;
  strLen=sprintf(&simHead[strPos],
                 "# Command line: %s",inp_argv[0]);
  strPos+=strLen;

  for(j=1; j<inp_argc; j++) {
    if(j < inp_argc-1) {
      strLen=sprintf(&simHead[strPos],
                     " %s",
                     inp_argv[j]);
    }
    else {
      strLen=sprintf(&simHead[strPos],
                     " %s\n%s\n",
                     inp_argv[j],
                     "# ---------------------------------------------------------------------------");
    }
    strPos+=strLen;
  }

  // 2. construct the mechanism information  
  strLen=sprintf(&simHead[strPos],
                 "%s%s%s\n%s%s\n%s%s\n%s%d\n%s%d\n%s%d\n%s\n",
                 "# Mechanism constants:\n",
                 "#   mechanism definition      [path] : ",
                 idt_ctrl->getMechFileName(),
                 "#   thermodynamics definition [path] : ",
                 idt_ctrl->getThermFileName(), 
                 "#   parser log definition     [path] : ",
                 idt_ctrl->getLogFileName(), 
                 "#   number of species            [#] : ",
                 nSpc, 
                 "#   number of reactions          [#] : ",
                 cvp->mechp->getNumReactions(),
                 "#   number of 1-way steps        [#] : ",
                 cvp->mechp->getNumSteps(),
		 //                 "#   max non-zeros in species J   [#] : ",
		 // nSpcTermsJac,
                 //"#   max non-zeros in const vol J [#] : ",
                 //cvp->sparseMtx->nNonZero,
                 "# ---------------------------------------------------------------------------");
  strPos+=strLen;
  // 3. construct the composition information
  strLen=sprintf(&simHead[strPos],
                 "%s\n",
                 "# Fuel composition (moles of species)/(moles of fuel):");
  strPos+=strLen;
  for(j=0; j<idt_ctrl->getNumFuelSpecies(); j++) {
    spcId = idt_ctrl->getFuelSpeciesId(j);
    strLen=sprintf(&simHead[strPos],
                   "#   %16s: %.18g\n",
                   cvp->mechp->getSpeciesName(spcId),
                   idt_ctrl->getFuelMoleFrac(spcId));
    strPos+=strLen;
  }
  strLen=sprintf(&simHead[strPos],
              "%s\n",
              "# Oxid(izer) composition (moles of species)/(moles of oxid):");
  strPos+=strLen;
  for(j=0; j<idt_ctrl->getNumOxidSpecies(); j++) {
    spcId = idt_ctrl->getOxidSpeciesId(j);
    strLen=sprintf(&simHead[strPos],
                   "#   %16s: %.18g\n",
                   cvp->mechp->getSpeciesName(spcId),
                   idt_ctrl->getOxidMoleFrac(spcId));
    strPos+=strLen;
  }
                   
  strLen=sprintf(&simHead[strPos],
                 "%s\n%s%.18g\n%s%.18g\n%s%.18g\n%s\n",
                 "# Global stoichiometry constants:",
                 "#   fuel (excess mols of atomic-O)/(mols of fuel) : ",
                 idt_ctrl->getFuelOxyBalance(),
                 "#   oxid (excess mols of atomic-O)/(mols of oxid) : ",
                 idt_ctrl->getOxidOxyBalance(),
                 "#   (mols of oxid)/(mols of fuel) at phi=1        : ",
                 idt_ctrl->getMoleOxidStoicPerFuel(),
                 "# ---------------------------------------------------------------------------");
  strPos+=strLen;


  // 4. construct the simulation information
  strLen=sprintf(&simHead[strPos],
                 "%s%s%.18g\n%s%.18g\n%s%s\n%s%.18g\n%s%.18g\n%s%d\n%s\n",
                 "# Simulation constants:\n",
                 "#   ode print time               [s] : ",
                 idt_ctrl->getPrintTime(),
                 "#   ode max time                 [s] : ",
                 idt_ctrl->getStopTime(),
                 "#   continue after IDT found   [y/n] : ",
                 (idt_ctrl->continueAfterIDT() != 0) ? "yes" : "no",
                 "#   reference temperature        [K] : ",
                 idt_ctrl->getRefTemp(),
                 "#   ignition delay delta Temp    [K] : ",
                 idt_ctrl->getDeltaTign(),
                 "#   total number of runs         [#] : ",
                 idt_ctrl->getRunTotal(),
                 "# ---------------------------------------------------------------------------");
  strPos+=strLen;

  // 5. construct the solver information
  strLen=sprintf(&simHead[strPos],
                 "%s%s%.18g\n%s%.18g\n%s%.18g\n%s%d\n%s%.18g\n%s%.18g\n%s%d\n%s%d\n%s%.18g\n%s%d\n%s%s\n%s%s\n%s%.18g\n%s%s\n%s%s\n%s%s\n%s\n",
                 "# Solver constants:\n",   
                 "#   ode relative tolerance            [-] : ",
                 idt_ctrl->getRTol(),
                 "#   ode absolute tolerance            [-] : ",
                 idt_ctrl->getATol(),
                 "#   ode max internal timestep         [s] : ",
                 idt_ctrl->getMaxInternalDt(),
                 "#   ode max internal steps            [#] : ",
                 idt_ctrl->getMaxInternalSteps(),
                 "#   linear solver rel error           [-] : ",
                 idt_ctrl->getEpsLin(),
                 "#   non-linear solver rel error       [-] : ",
                 idt_ctrl->getNlConvCoeff(),
                 "#   cvode primary failures allowed    [#] : ",
                 idt_ctrl->getMaxPrimaryCvodeFails(),
                 "#   cvode secondary failures allowed  [#] : ",
                 idt_ctrl->getMaxSecondaryCvodeFails(),
                 "#   secondary failure P threshold     [-] : ",
                 idt_ctrl->getSafetyThreshold(),
                 "#   max Krylov dimension              [#] : ",
                 idt_ctrl->getKrylovDim(),
                 "#   use fake P update               [y/n] : ",
                 (idt_ctrl->getUpdate() ? "yes" : "no"),
                 "#   use incomplete LU factorization [y/n] : ",
                 (idt_ctrl->getILU() ? "yes" : "no"),
                 "#   partial pivot threshold           [-] : ",
                 idt_ctrl->getPartialPivotThresh(),
                 "#   P threshold type code             [#] : ",
                 precondThresholdName[idt_ctrl->getThreshType()-1],
                 "#   permutation type code             [#] : ",
                 permutationTypeName[idt_ctrl->getPermutationType()-1],
                 "#   exp() vector                      [#] : ",
                 zerork::expType,
                 "# ---------------------------------------------------------------------------");
  if(strPos >= maxLen) {
    printf("ERROR: the simulation header info (%d bytes)\n",strPos);
    printf("       is larger than the buffer size (%d bytes).\n",maxLen);
    printf("       Exiting Now!\n");
    exit(-1);
  }
  
  return strPos; 
}

int getIdtLine(const double initMoleFrac[],
               const idt_sweep_params *idt_ctrl,
               const double idt,
               const double runTime,
               const int nCvodeFails,
               const int maxBuffer,
               char *idtLine)
{
  int j;
  int nTrack = idt_ctrl->getNumTrackSpecies();
  int strPos,strLen;

  strPos=0;
  strLen=sprintf(&idtLine[strPos],
                 "%8d  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14d",
                 idt_ctrl->getRunId(),
		 idt_ctrl->getInitTemp(),
                 idt_ctrl->getInitPres(),
                 idt_ctrl->getPhi(),
                 idt_ctrl->getEgr(),
                 idt_ctrl->getThresh(),
                 idt_ctrl->getDensity(),
                 idt,
                 runTime,
                 nCvodeFails);
  strPos+=strLen;
  
  for(j=0; j<nTrack; j++) {
    strLen=sprintf(&idtLine[strPos],
                   "   %14.7e",
                   initMoleFrac[idt_ctrl->getTrackSpeciesId(j)]);
    strPos+=strLen;
  }
  strLen=sprintf(&idtLine[strPos],"\n");
  strPos+=strLen;

  if(strPos >= maxBuffer) {
    printf("ERROR: the IDT file data line (%d bytes)\n",strPos);
    printf("       is larger than the buffer size (%d bytes).\n",maxBuffer);
    printf("       Exiting Now!\n");
    exit(-1);
  }

  return strPos;
}
int getIdtHeader(const idt_sweep_params *idt_ctrl,
                 const ODEParams_CMT *cvp,
                 const int maxBuffer,
                 char *idtHeader)
{
  int j;
  int nTrack = idt_ctrl->getNumTrackSpecies();
  int strPos,strLen;

  strPos=0;
  strLen=sprintf(&idtHeader[strPos],
                 "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
                 "# Column  1: [#]      run id",
                 "# Column  2: [K]      initial temperature",
                 "# Column  3: [Pa]     initial pressure",
                 "# Column  4: [-]      equivalence ratio",
                 "# Column  5: [-]      EGR fraction",
                 "# Column  6: [-]      preconditioner threshold",
                 "# Column  7: [kg/m^3] initial density",
                 "# Column  8: [s]      ignition delay time (defined delta T)",
                 "# Column  9: [s]      execution wall clock time",
                 "# Column 10: [#]      number of cvode errors");
  strPos+=strLen;

  for(j=0; j<nTrack; j++) {
    strLen=sprintf(&idtHeader[strPos],
                   "# Column %2d: [-]      initial mole frac of %s\n",
                   j+11,
                   cvp->mechp->getSpeciesName(idt_ctrl->getTrackSpeciesId(j)));
    strPos+=strLen;
  }
  strLen=sprintf(&idtHeader[strPos],"# ----------------------------------------------------------------------------\n");
  strPos+=strLen;  

  strLen=sprintf(&idtHeader[strPos],
                 "%7s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |%14s |",
                 "# run id",
                 "Temp [K]",
                 "Pres [Pa]",
                 "phi [-]",
                 "EGR frac [-]",
                 "P thresh [-]",
                 "dens [kg/m3]",
                 "IDT [s]",
                 "sim time [s]",
                 "cvode err [#]");
  strPos+=strLen;

  for(j=0; j<nTrack; j++) {
    strLen=sprintf(&idtHeader[strPos],
                   " mlfrc %8s |",
                   cvp->mechp->getSpeciesName(idt_ctrl->getTrackSpeciesId(j)));
    strPos+=strLen;
  }
  strLen=sprintf(&idtHeader[strPos],"\n");
  strPos+=strLen;  

  if(strPos >= maxBuffer) {
    printf("ERROR: the IDT file header info (%d bytes)\n",strPos);
    printf("       is larger than the buffer size (%d bytes).\n",maxBuffer);
    printf("       Exiting Now!\n");
    exit(-1);
  }

  return strPos;
}
 


/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
