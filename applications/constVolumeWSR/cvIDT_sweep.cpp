#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const int MAX_SPECNAME_LEN=256;
const int MAX_FLOAT_LEN=256;
const int MAX_THIST_LINE=20000;

#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros

#ifdef SUNDIALS2
#include <cvode/cvode_spgmr.h>      // prototypes & constants for CVSPGMR
#elif SUNDIALS3
#include <cvode/cvode_spils.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#elif SUNDIALS4
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif

#include <CKconverter/CKReader.h>
#include <zerork/constants.h>
#include <zerork/constants_api.h>

#include "matrix_funcs.h"
#include "cv_param_sparse.h"
#include "ode_funcs.h"
#ifdef SPIFY
#include "sweep_util_yml.h"
#else
#include "sweep_util_formatted.h"
#endif
#include "utility_funcs.h"


int setMoleFracFromFile(zerork::mechanism &mechInp, const char *fileName,
			double moleFrac[]);

static int check_flag(void *flagvalue, const char *funcname, int opt);

void getTimeHistLine_full(const double currTime,
                          const double *sysState,
                          const idt_sweep_params *idt_ctrl,
                          const cv_param *cvp,
                          const double wallTime,
                          char *thistLine);

void getTimeHistHeader_full(const idt_sweep_params *idt_ctrl,
                            const cv_param *cvp,
                            char *thistHead);
int getIdtLine(const double initMoleFrac[],
               const idt_sweep_params *idt_ctrl,
               const double idt,
               const double runTime,
               const int nCvodeFails,
               const int maxBuffer,
               char *idtLine);
int getIdtHeader(const idt_sweep_params *idt_ctrl,
                 const cv_param *cvp,
                 const int maxBuffer,
                 char *idtHeader);


int getSimHeader_full(const int inp_argc,
                      char **inp_argv,
                      const idt_sweep_params *idt_ctrl,
                      const cv_param *cvp,
                      const int max,
                      char *simHead);


void cvReactor(int inp_argc, char **inp_argv)
{
  double idt,tmax,dtprint;
  int isIdtFound;
  double nRunsPerThresh;

  FILE *idtFilePtr,*thistFilePtr;
  int j,k;
  int nTrackSpec,nSpc,nState,nStep;
  double tnext,tcurr;
  // cvode variables
  N_Vector systemState;
  cv_param systemParam;
  void *cvode_mem;
#ifdef SUNDIALS3
  SUNLinearSolver LS;
#elif SUNDIALS4
  SUNLinearSolver LS;
  SUNNonlinearSolver NLS;
#endif
  int flag;

  char thistLine[MAX_THIST_LINE];
  int did_cvodeFail;
  int isBadStep;

  // timing data
  double startTime,stopTime, simTime;
//  double setupZeroRK,setupCVode,setupJterm;

  if(inp_argc != 2)
    {
      printf("ERROR: incorrect command line usage\n");
      printf("       use instead %s <idt sweep input>\n",inp_argv[0]);
      fflush(stdout); exit(-1);
    }

  startTime=getHighResolutionTime();
  // sweep controls
  idt_sweep_params idt_ctrl(inp_argv[1]);

  stopTime=getHighResolutionTime();
//  setupZeroRK=stopTime-startTime;

  nSpc=idt_ctrl.getNumSpecies();
  nState=idt_ctrl.getNumSpecies()+1;
  nStep=idt_ctrl.getNumSteps();

  std::vector<double> avgTime(idt_ctrl.getNumThreshRuns(), 0.0);
  std::vector<double> minTime(idt_ctrl.getNumThreshRuns(), 1.0e300);
  std::vector<double> maxTime(idt_ctrl.getNumThreshRuns(), -1.0e300);

  systemState = N_VNew_Serial(nState);
  double* massFracPtr=NV_DATA_S(systemState); // caution: assumes realtype == double
  std::vector<double> moleFracCurr(nSpc);
  std::vector<double> moleFracInit(nSpc);

  // set up the system parameters
  systemParam.nSpc=idt_ctrl.getNumSpecies();
  systemParam.minMassFrac=1.0e-30;
  systemParam.sqrtUnitRnd=sqrt(UNIT_ROUNDOFF);
  systemParam.Tref=idt_ctrl.getRefTemp();
  systemParam.mech=idt_ctrl.getMechPtr();

  // don't forget to reset densities for each calculations
  //systemParam.Dens=gasMech->getDensityFromTPY(tempSweep[0],pres0,massFracPtr);
  //systemParam.invDens=1.0/systemParam.Dens;
  systemParam.netProd    =(double *)malloc(sizeof(double)*nSpc);
  systemParam.Energy     =(double *)malloc(sizeof(double)*nSpc);
  systemParam.CvMass     =(double *)malloc(sizeof(double)*nSpc);
  systemParam.molWt      =(double *)malloc(sizeof(double)*nSpc);
  systemParam.invMolWt   =(double *)malloc(sizeof(double)*nSpc);
  systemParam.fwdROP     =(double *)malloc(sizeof(double)*nStep);
  systemParam.createRate =(double *)malloc(sizeof(double)*nSpc);
  systemParam.destroyRate=(double *)malloc(sizeof(double)*nSpc);
  systemParam.conc       =(double *)malloc(sizeof(double)*nSpc);

  // set constant parameters
  systemParam.mech->getMolWtSpc(systemParam.molWt);
  for(j=0; j<nSpc; j++)
    {systemParam.invMolWt[j]=1.0/systemParam.molWt[j];}

  // set up sparse jacobian matrix class
  startTime=getHighResolutionTime();
  systemParam.sparseMtx=
        (Jsparse *)alloc_Jsparse(*systemParam.mech,0.0,idt_ctrl.getILU(),
                                 idt_ctrl.getUpdate(),idt_ctrl.getThreshType(),
                                 idt_ctrl.getPartialPivotThresh(), idt_ctrl.getPermutationType());
  stopTime=getHighResolutionTime();
//  setupJterm=stopTime-startTime;

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration */
  startTime=getHighResolutionTime();
#ifdef SUNDIALS4
  cvode_mem = CVodeCreate(CV_BDF);
#else
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
#endif
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) exit(-1);
  systemParam.cvodeMemPtr=cvode_mem;

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag=CVodeInit(cvode_mem, const_vol_wsr, 0.0, systemState);
  if (check_flag(&flag, "CVodeInit", 1)) exit(-1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  flag = CVodeSStolerances(cvode_mem,idt_ctrl.getRTol(),idt_ctrl.getATol());
  if (check_flag(&flag, "CVodeSStolerances", 1)) exit(-1);

  /* Call CVodeRootInit to specify the root function with 1 component */
  flag = CVodeRootInit(cvode_mem, 1, tempRootFunc);
  //flag = CVodeRootInit(cvode_mem, 2, tempRoot2);
  if (check_flag(&flag, "CVodeRootInit", 1)) exit(-1);

  /* Set the pointer to user-defined data */
  flag = CVodeSetUserData(cvode_mem, &systemParam);
  if(check_flag(&flag, "CVodeSetUserData", 1)) exit(-1);

  /* Call CVSpgmr to specify the linear solver CVSPGMR
     with left preconditioning and the maximum Krylov dimension maxl */
#ifdef SUNDIALS2
  flag = CVSpgmr(cvode_mem, PREC_LEFT, idt_ctrl.getKrylovDim());
  if(check_flag(&flag, "CVSpgmr", 1)) exit(-1);

  flag = CVSpilsSetGSType(cvode_mem, MODIFIED_GS);
  if(check_flag(&flag, "CVSpilsSetGSType", 1)) exit(-1);

  /* Set preconditioner setup and solve routines Precond and PSolve,
     and the pointer to the user-defined block data */
  flag = CVSpilsSetPreconditioner(cvode_mem, jac_full_prec_setup,
  				  jac_full_prec_solveV3);
  if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) exit(-1);

#elif SUNDIALS3
  LS = SUNSPGMR(systemState, PREC_LEFT, idt_ctrl.getKrylovDim());

  flag = CVSpilsSetLinearSolver(cvode_mem, LS);
  if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) exit(-1);

  flag = SUNSPGMRSetGSType(LS, MODIFIED_GS);
  if(check_flag(&flag, "SUNSPGMRSetGSType", 1)) exit(-1);
    /* Set preconditioner setup and solve routines Precond and PSolve,
     and the pointer to the user-defined block data */
  flag = CVSpilsSetPreconditioner(cvode_mem, jac_full_prec_setup,
  				  jac_full_prec_solveV3);
  if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) exit(-1);

#elif SUNDIALS4
  NLS = SUNNonlinSol_Newton(systemState);
  flag = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if(check_flag(&flag, "CVodeSetNonlinearSolver", 1)) exit(-1);
  LS = SUNLinSol_SPGMR(systemState, PREC_LEFT, idt_ctrl.getKrylovDim());
  flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
  if(check_flag(&flag, "CVodeSetLinearSolver", 1)) exit(-1);

  /* Set preconditioner setup and solve routines Precond and PSolve,
     and the pointer to the user-defined block data */
  flag = CVodeSetPreconditioner(cvode_mem, jac_full_prec_setup,
                                jac_full_prec_solveV3);
  if(check_flag(&flag, "CVodeSetPreconditioner", 1)) exit(-1);
#endif

//  flag = CVSpilsSetJacTimesVecFn(cvode_mem, sparse_jac_v);
//  if(check_flag(&flag, "CVSpilsSetJacTimesVecFn", 1)) exit(-1);

  /* Set the maximum number of internal steps per CVode call and the maximum
   * allowable internal steps. */
  flag = CVodeSetMaxNumSteps(cvode_mem, idt_ctrl.getMaxInternalSteps());
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) exit(-1);
  flag = CVodeSetMaxStep(cvode_mem, idt_ctrl.getMaxInternalDt());
  if (check_flag(&flag, "CVodeSetMaxStep", 1)) exit(-1);

  flag = CVodeSetNonlinConvCoef(cvode_mem, idt_ctrl.getNlConvCoeff()); // Default [0.1]
#if defined SUNDIALS2 || defined SUNDIALS3
  flag = CVSpilsSetEpsLin(cvode_mem, idt_ctrl.getEpsLin());    // Default [0.05]
#else
  flag = CVodeSetEpsLin(cvode_mem, idt_ctrl.getEpsLin());    // Default [0.05]
#endif
  stopTime=getHighResolutionTime();
//  setupCVode=stopTime-startTime;
  nTrackSpec = idt_ctrl.getNumTrackSpecies();

  // ready for integration

  // open ignition delay and time history files for write
  idtFilePtr=fopen(idt_ctrl.getIdtFileName(),"w");
  if(idtFilePtr==NULL)
    {
      printf("ERROR: could not open output file %s for write\n",
	     idt_ctrl.getIdtFileName());
      exit(-1);
    }
  thistFilePtr=fopen(idt_ctrl.getTHistFileName(),"w");
  if(thistFilePtr==NULL)
    {
      printf("ERROR: could not open output file %s for write\n",
	     idt_ctrl.getTHistFileName());
      exit(-1);
    }
  // create the simulation header stored in a character array
  getSimHeader_full(inp_argc,
                    &inp_argv[0],
                    &idt_ctrl,
                    &systemParam,
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
  getTimeHistHeader_full(&idt_ctrl, &systemParam, thistLine);
  // write the time history header to stdout and the time history file
  printf("%s",thistLine);
  fflush(stdout);
  fprintf(thistFilePtr,"%s",thistLine);
  fflush(thistFilePtr);

  // create the idt file header stored in a character array
  getIdtHeader(&idt_ctrl, &systemParam, MAX_THIST_LINE, thistLine);
  // write the idt header to the idt file
  fprintf(idtFilePtr,"%s",thistLine);
  fflush(idtFilePtr);

  tmax    = idt_ctrl.getStopTime();
  dtprint = idt_ctrl.getPrintTime();

  for(j=0; j<idt_ctrl.getRunTotal(); j++)
    {
      // initialize the run counters
      systemParam.sparseMtx->reduceNNZ=0; // won't be set until the first J
      systemParam.sparseMtx->LUnnz = 0; // won't be set until the first J
      systemParam.sparseMtx->fillFactor= 100.; // won't be set until the first J
      systemParam.prevNumErrTestFails = 0;
      systemParam.nFunc=0;
      systemParam.nJacSetup=0;
      systemParam.nJacFactor=0;
      systemParam.nBackSolve=0;
      systemParam.nJacRescale=0;
      systemParam.nColPerm=0;
      systemParam.colPermTime   = 0;
      systemParam.jacFactorTime = 0;
      systemParam.backsolveTime = 0;
      systemParam.jacSetupTime = 0;
      systemParam.funcTime = 0;
      did_cvodeFail=0;
      isIdtFound=0;

      // reset the ignition delay threshold
      systemParam.tempRoot=idt_ctrl.getDeltaTign()+idt_ctrl.getInitTemp();

      // reset the preconditioner threshold
      change_JsparseThresh(systemParam.sparseMtx,idt_ctrl.getThresh());

      // reset the mass fraction density and relative volume
      idt_ctrl.getInitMassFrac(massFracPtr);
      systemParam.Dens=idt_ctrl.getDensity();
      systemParam.invDens=1.0/systemParam.Dens;
      NV_Ith_S(systemState,nSpc)=
	idt_ctrl.getInitTemp()/idt_ctrl.getRefTemp();

      // reset the time
      tcurr=0.0;
      tnext=dtprint;
      idt=tmax;

      // reinitialize cvode
      flag = CVodeReInit(cvode_mem, tcurr, systemState);
      isBadStep = 0;

      getTimeHistLine_full(tcurr,NV_DATA_S(systemState),&idt_ctrl,
                           &systemParam,0.0,thistLine);
      printf("%s",thistLine);
      fflush(stdout);
      fprintf(thistFilePtr,"%s",thistLine);
      fflush(thistFilePtr);
      startTime=getHighResolutionTime();

      while(tcurr<tmax)
  	{
          if(idt_ctrl.oneStep()) {
            flag = CVode(cvode_mem, tnext, systemState, &tcurr, CV_ONE_STEP);
          }
          else {
  	    flag = CVode(cvode_mem, tnext, systemState, &tcurr, CV_NORMAL);
          }

  	  if (check_flag(&flag, "CVODE ERROR", 1)) {

	    isBadStep = 1;
	    did_cvodeFail++;
            if(did_cvodeFail <= idt_ctrl.getMaxPrimaryCvodeFails()) {

	      printf("WARNING: attempting CVodeReInit at t=%.18g [s]\n",tcurr);
	      printf("         after cvode failure count = %d\n",
              did_cvodeFail);

	      flag = CVodeReInit(cvode_mem, tcurr, systemState);
            }
            else if(did_cvodeFail <= idt_ctrl.getMaxPrimaryCvodeFails() +
                    idt_ctrl.getMaxSecondaryCvodeFails()) {
	      printf("WARNING: resetting the preconditioner threshold to %g\n",
                     idt_ctrl.getSafetyThreshold());
	      // reset the preconditioner threshold
	      change_JsparseThresh(systemParam.sparseMtx,
                                   idt_ctrl.getSafetyThreshold());

	      printf("WARNING: attempting CVodeReInit at t=%.18g [s]\n",tcurr);
	      printf("         after cvode failure count = %d\n",
                     did_cvodeFail);

              flag = CVodeReInit(cvode_mem, tcurr, systemState);
	    }
            else {
              printf("ERROR: idt calculation failed to recover after\n");
	      printf("       %d re-initializations.  Halting Now!\n",
                     did_cvodeFail-1);
		exit(-1);
            }
	  }

          getTimeHistLine_full(tcurr,
                               NV_DATA_S(systemState),
                               &idt_ctrl,
                               &systemParam,
                               getHighResolutionTime()-startTime,
                               thistLine);
          printf("%s",thistLine);
          fflush(stdout);
          fprintf(thistFilePtr,"%s",thistLine);
          fflush(thistFilePtr);

  	  if(flag == CV_ROOT_RETURN) {
            idt=tcurr;
            isIdtFound=1;
            if(idt_ctrl.continueAfterIDT() == 0) {
              tcurr=2.0*tmax; // advance time past the stopping criteria
	    }
          }
  	  else if(!idt_ctrl.oneStep() && isBadStep==0)
  	    {tnext+=dtprint;}
          else if(tcurr >= tnext) {
            tnext += dtprint;
            flag = CVodeReInit(cvode_mem, tcurr, systemState);
          }
	}

      if(idt_ctrl.dumpJacobian()) {
          char jacFileName[32];
          snprintf(jacFileName,32,"jacobian%03d.txt",idt_ctrl.getRunId());
          FILE* jacFile=fopen(jacFileName,"w");
          print_sp_matrix(jacFile,systemParam.nSpc+1,systemParam.nSpc+1,
                          systemParam.sparseMtx->mtxColSum,
                          systemParam.sparseMtx->mtxRowIdx,
                          systemParam.sparseMtx->mtxData);
          fclose(jacFile);
      }

      stopTime=getHighResolutionTime();
      simTime=stopTime-startTime;
      // add new-line to separate the time histories
      printf("\n");
      fprintf(thistFilePtr,"\n");

      idt_ctrl.getInitMoleFrac(&moleFracInit[0]);
      getIdtLine(&moleFracInit[0],&idt_ctrl,idt,simTime,did_cvodeFail,
                 MAX_THIST_LINE,thistLine);
      fprintf(idtFilePtr,"%s",thistLine);
      fflush(idtFilePtr);

      avgTime[idt_ctrl.getThreshId()]+=simTime;
      if(simTime > maxTime[idt_ctrl.getThreshId()])
	{maxTime[idt_ctrl.getThreshId()]=simTime;}
      if(simTime < minTime[idt_ctrl.getThreshId()])
	{minTime[idt_ctrl.getThreshId()]=simTime;}

      idt_ctrl.incrementRunId();
    }

  // report the threshold statistics
  fprintf(idtFilePtr,"\n\n");
  fprintf(idtFilePtr,"# -----------------------------------------------------------------------------\n");
  fprintf(idtFilePtr,"# Column 1: [-] preconditioner threshold\n");
  fprintf(idtFilePtr,"# Column 2: [s] average simulation time\n");
  fprintf(idtFilePtr,"# Column 3: [s] minimum simulation time\n");
  fprintf(idtFilePtr,"# Column 4: [s] maximum simulation time\n");
  nRunsPerThresh=(double)(idt_ctrl.getRunTotal()/
			  idt_ctrl.getNumThreshRuns());
  for(j=0; j<idt_ctrl.getNumThreshRuns(); j++)
    {
      avgTime[j]/=nRunsPerThresh;
      fprintf(idtFilePtr,"%14.7e  %14.7e  %14.7e  %14.7e\n",
	      idt_ctrl.getThreshFromId(j),avgTime[j],minTime[j],maxTime[j]);
    }
  fflush(idtFilePtr);

  // clean up procedures
  fclose(idtFilePtr);
  fclose(thistFilePtr);
 // Free ode system parameter memory
  free(systemParam.molWt);
  free(systemParam.invMolWt);
  free(systemParam.netProd);
  free(systemParam.fwdROP);
  free(systemParam.createRate);
  free(systemParam.destroyRate);
  free(systemParam.conc);
  free(systemParam.Energy);
  free(systemParam.CvMass);

  free_Jsparse(systemParam.sparseMtx);
  N_VDestroy_Serial(systemState);

  // Free integrator memory
  CVodeFree(&cvode_mem);
#if defined SUNDIALS3 || defined SUNDIALS4
  SUNLinSolFree(LS);
#endif
#if defined SUNDIALS4
  SUNNonlinSolFree(NLS);
#endif

//  printf("Begin.......delete of idt_sweep_params class\n"); fflush(stdout);
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
void getTimeHistLine_full(const double currTime,
                          const double *sysState,
                          const idt_sweep_params *idt_ctrl,
                          const cv_param *cvp,
                          const double wallTime,
                          char *thistLine)
{
  int j;
  int nSpc   = cvp->mech->getNumSpecies();
  int nTrack = idt_ctrl->getNumTrackSpecies();
  int strPos,strLen;
  double currPres,currTemp, initTemp;
  double *currMoleFrac;
  long int nsteps,netfails,nniters,nncfails;

  double molWt, cvMixPerMass, intEnergyMixPerMass, intEnergyMixPerMassInit;

  nsteps=netfails=nniters=nncfails=0;

  CVodeGetNumSteps(cvp->cvodeMemPtr, &nsteps);
  CVodeGetNumErrTestFails(cvp->cvodeMemPtr, &netfails);
  CVodeGetNumNonlinSolvIters(cvp->cvodeMemPtr, &nniters);
  CVodeGetNumNonlinSolvConvFails(cvp->cvodeMemPtr, &nncfails);

  currMoleFrac = new double[nSpc];

  initTemp = idt_ctrl->getInitTemp();
  currTemp = sysState[nSpc]*idt_ctrl->getRefTemp();
  currPres = cvp->mech->getPressureFromTVY(currTemp,
                                           cvp->invDens,
                                           &sysState[0]);

  molWt=cvp->mech->getMolWtMixFromY(&sysState[0]);
  // currMoleFrac[] used as temporary storage for species specific heat
  cvMixPerMass=cvp->mech->getMassCvFromTY(currTemp,
                                          &sysState[0],
                                          &currMoleFrac[0]);
  // currMoleFrac[] used as temporary storage for species internal energy
  intEnergyMixPerMass=cvp->mech->getMassIntEnergyFromTY(currTemp,
                                                        &sysState[0],
                                                        &currMoleFrac[0]);
  // currMoleFrac[] used as temporary storage for species internal energy
  intEnergyMixPerMassInit=cvp->mech->getMassIntEnergyFromTY(initTemp,
                                                            &sysState[0],
                                                            &currMoleFrac[0]);
  // currMoleFrac[] is recalculated with mole fraction
  cvp->mech->getXfromY(&sysState[0],currMoleFrac);

  // sprintf return the number of characters written excluding the '\0'
  // string termination character
  if(idt_ctrl->longOutput()) {
    strLen=sprintf(thistLine,
  		 "%8d  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e  %24.16e",
		 idt_ctrl->getRunId(),
		 currTime,
		 currTemp,
		 currPres,
		 1.0/cvp->invDens,
                 molWt,
                 cvMixPerMass,
                 intEnergyMixPerMass-intEnergyMixPerMassInit);
    strPos=strLen;
    for(j=0; j<nTrack; j++)
      {
        strLen=sprintf(&thistLine[strPos],
                "  %24.16e",
                currMoleFrac[idt_ctrl->getTrackSpeciesId(j)]);
        strPos+=strLen;
      }
  } else {
    strLen=sprintf(thistLine,
  		 "%8d  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
		 idt_ctrl->getRunId(),
		 currTime,
		 currTemp,
		 currPres,
		 1.0/cvp->invDens,
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
  }

  strLen=sprintf(&thistLine[strPos],
                 "  %6ld  %6ld  %6ld  %6ld  %6d  %6d  %6d  %6d  %6d  %6d  %6d  %6d",
                 nsteps,
                 netfails,
                 nniters,
                 nncfails,
                 cvp->nJacSetup,
                 cvp->nColPerm,
                 cvp->nJacFactor,
                 cvp->nBackSolve,
                 cvp->nFunc,
	         cvp->sparseMtx->nNonZero,
	         cvp->sparseMtx->reduceNNZ,
	         cvp->sparseMtx->LUnnz);
  strPos+=strLen;
  strLen=sprintf(&thistLine[strPos],
                 "  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e\n",
                 cvp->jacSetupTime,
                 cvp->colPermTime,
                 cvp->jacFactorTime,
                 cvp->backsolveTime,
                 cvp->funcTime,wallTime);

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
                            const cv_param *cvp,
                            char *thistHead)
{
  int j;
  int nSpc   = cvp->mech->getNumSpecies();
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
                   cvp->mech->getSpeciesName(idt_ctrl->getTrackSpeciesId(j)));
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
                 "# Column %2d: [#]      calls to the precond. setup\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [#]      calls to the column permutation\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [#]      calls to the precond. factorization\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [#]      calls to the precond. back solve\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [#]      calls to the RHS function\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [#]      max non-zero terms in precond.\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [#]      current non-zero terms in precond.\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [#]      current non-zero terms in LU factors\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [s]      wall clock time - precond. setup\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [s]      wall clock time - column permut.\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [s]      wall clock time - precond. factor.\n",
                 colNum++);
  strPos+=strLen;
  strLen=sprintf(&thistHead[strPos],
                 "# Column %2d: [s]      wall clock time - precond. back solve\n",
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
  if(idt_ctrl->longOutput()) {
    strLen=sprintf(&thistHead[strPos], "# run id  %24s  %24s  %24s  %24s  %24s  %24s  %24s",
                            "sim time","temp", "press", "dens","mol wt", "cv", "E release");
    strPos+=strLen;
    for(j=0; j<nTrack; j++) {
      strLen=sprintf(&thistHead[strPos],
                     " mlfrc%20s",
                     cvp->mech->getSpeciesName(idt_ctrl->getTrackSpeciesId(j)));
      strPos+=strLen;
    }
  } else {
    strLen=sprintf(&thistHead[strPos],
                   "# run id        sim time            temp           press            dens          mol wt              cv       E release");
    strPos+=strLen;
    for(j=0; j<nTrack; j++) {
      strLen=sprintf(&thistHead[strPos],
                     " mlfrc%10s",
                     cvp->mech->getSpeciesName(idt_ctrl->getTrackSpeciesId(j)));
      strPos+=strLen;
    }
  }
  strLen=sprintf(&thistHead[strPos],
                 "   nstep   nertf  nnlsit  nnlcvf prsetup ncolprm nfactor backslv  intfcl     nnz    rnnz   lunnz    jac setup tm     col perm tm      jac fac tm    backsolve tm     rhs func tm   simulation tm\n");
  strPos+=strLen;
}

// void getSimHeader_full(...)
//
//
int getSimHeader_full(const int inp_argc,
                      char **inp_argv,
                      const idt_sweep_params *idt_ctrl,
                      const cv_param *cvp,
                      const int maxLen,
                      char *simHead)
{
  int j,spcId;
  int nSpc   = cvp->mech->getNumSpecies();
  int nSpcTermsJac = cvp->sparseMtx->nNonZero-2*nSpc+1;
  int nTrack = idt_ctrl->getNumTrackSpecies();
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

  {
    zerork::ZeroRKConstants zrk_constants;
    // 1.5 construct git commit information
    strLen=sprintf(&simHead[strPos],
                   "%s%s%s\n%s%s\n%s%s\n%s%s\n%s\n",
                   "# Zero-RK git commit information:\n",
                   "#   commit id                 [hash] : ",
                   zrk_constants.GetCommitId(),
                   "#   commit time-stamp         [date] : ",
                   zrk_constants.GetCommitTimestamp(),
                   "#   branch name               [word] : ",
                   zrk_constants.GetBranchName(),
                   "#   exponential type          [word] : ",
                   zrk_constants.GetExpType(),
                   "# ---------------------------------------------------------------------------");
    strPos+=strLen;
  }


  // 2. construct the mechanism information
  strLen=sprintf(&simHead[strPos],
                 "%s%s%s\n%s%s\n%s%s\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n%s\n",
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
                 cvp->mech->getNumReactions(),
                 "#   number of 1-way steps        [#] : ",
                 cvp->mech->getNumSteps(),
                 "#   max non-zeros in species J   [#] : ",
                 nSpcTermsJac,
                 "#   max non-zeros in const vol J [#] : ",
                 cvp->sparseMtx->nNonZero,
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
                   cvp->mech->getSpeciesName(spcId),
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
                   cvp->mech->getSpeciesName(spcId),
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
                 "%s%s%.18g\n%s%.18g\n%s%.18g\n%s%d\n%s%.18g\n%s%.18g\n%s%d\n%s%d\n%s%.18g\n%s%d\n%s%s\n%s%s\n%s%.18g\n%s%s\n%s%s\n%s\n",
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
                 "%8d  %23.16e  %23.16e  %23.16e  %23.16e  %23.16e  %23.16e  %23.16e  %23.16e  %14d",
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
                   "   %23.16e",
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
                 const cv_param *cvp,
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
                   cvp->mech->getSpeciesName(idt_ctrl->getTrackSpeciesId(j)));
    strPos+=strLen;
  }
  strLen=sprintf(&idtHeader[strPos],"# ----------------------------------------------------------------------------\n");
  strPos+=strLen;

  strLen=sprintf(&idtHeader[strPos],
                 "%7s |%23s |%23s |%23s |%23s |%23s |%23s |%23s |%23s |%14s |",
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
                   " mlfrc %17s |",
                   cvp->mech->getSpeciesName(idt_ctrl->getTrackSpeciesId(j)));
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
