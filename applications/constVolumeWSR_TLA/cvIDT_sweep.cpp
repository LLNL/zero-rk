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
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <numeric>

#ifdef ZERORK_MPI
#include "mpi.h"
#endif

const int MAX_SPECNAME_LEN=256;
const int MAX_FLOAT_LEN=256;

//#include <culapack.h>

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

using zerork::getHighResolutionTime;

typedef struct
{
  int rxnId;
  double relSens;
} rxnSens_t;

int compare_rxnSens_t(const void *A, const void *B);

int setMoleFracFromFile(zerork::mechanism &mechInp, const char *fileName,
			double moleFrac[]);

static int check_flag(void *flagvalue, const char *funcname, int opt);
#ifdef ZERORK_MPI
static void mpi_sendrecv_string(int from, int to, std::string& in, std::string* out, MPI_Comm comm);
#endif

void getTimeHistLine_full(const double currTime,
                          const double *sysState,
                          const idt_sweep_params *idt_ctrl,
                          const cv_param *cvp,
                          const double wallTime,
                          std::string *thistLine);

void getTimeHistHeader_full(const idt_sweep_params *idt_ctrl,
                            const cv_param *cvp,
                            std::string *thistHead);
void getIdtLine(const double initMoleFrac[],
               const idt_sweep_params *idt_ctrl,
               const std::vector<double> idts,
               const double runTime,
               const int nCvodeFails,
               std::string *idtLine);
void getIdtHeader(const idt_sweep_params *idt_ctrl,
                 const cv_param *cvp,
                 std::string *idtHeader);


void getSimHeader_full(const int inp_argc,
                      char **inp_argv,
                      const idt_sweep_params *idt_ctrl,
                      const cv_param *cvp,
                      std::string *simHead);


void cvReactor(int inp_argc, char **inp_argv)
{
  int mpi_rank = 0;
  int mpi_size = 1;
#ifdef ZERORK_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif


  double tmax,dtprint;
  std::vector<std::vector<double>> idts;
  int isIdtFound;
  double nRunsPerThresh;

  FILE *idtFilePtr,*thistFilePtr;
  int i,j,k;
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

  std::string thistLine;
  int did_cvodeFail;
  int isBadStep;

  // timing data
  double startTime,stopTime, simTime;
//  double setupZeroRK,setupCVode,setupJterm;

  if(inp_argc != 2)
    {
      if(mpi_rank == 0)
      {
      printf("ERROR: incorrect command line usage\n");
      printf("       use instead %s <idt sweep input>\n",inp_argv[0]);
          fflush(stdout);
      }
#ifdef ZERORK_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      exit(-1);
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
#ifdef ZERORK_MPI
  std::vector<double> global_avgTime(idt_ctrl.getNumThreshRuns(), 0.0);
  std::vector<double> global_minTime(idt_ctrl.getNumThreshRuns(), 1.0e300);
  std::vector<double> global_maxTime(idt_ctrl.getNumThreshRuns(), -1.0e300);
#endif

  systemState = N_VNew_Serial(nState);
  double* massFracPtr=NV_DATA_S(systemState); // caution: assumes realtype == double
  std::vector<double> moleFracCurr(nSpc);
  std::vector<double> moleFracInit(nSpc);

  // temp vector to checkpoint state
  int nMaxTimeSteps = 10000; //?
  std::vector<double> systemStateCheckpoint, timeCheckpoint;
  systemStateCheckpoint.assign(nMaxTimeSteps*nState, 0.0);
  timeCheckpoint.assign(nMaxTimeSteps, 0.0);
  double dTdt;

  // set up the system parameters
  systemParam.nSpc=idt_ctrl.getNumSpecies();
  systemParam.minMassFrac=1.0e-30;
  systemParam.sqrtUnitRnd=sqrt(UNIT_ROUNDOFF);
  systemParam.Tref=idt_ctrl.getRefTemp();
  systemParam.mech=idt_ctrl.getMechPtr();
  systemParam.numTempRoots = idt_ctrl.getTemperatureDeltas().size();
  systemParam.tempRoots.resize(systemParam.numTempRoots);

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
  // previous timestep Jacobian
  systemParam.sparseMtxOld=
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

  /* Call CVodeRootInit to specify the root function */
  flag = CVodeRootInit(cvode_mem, systemParam.numTempRoots, tempRootFunc);
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

  if(mpi_rank == 0) {
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
                      &thistLine);
    // write the simulation  header to stdout, idt file and thist file
    printf("%s",thistLine.c_str());
    fflush(stdout);
    fprintf(thistFilePtr,"%s",thistLine.c_str());
    fflush(thistFilePtr);
    fprintf(idtFilePtr,"%s",thistLine.c_str());
    fflush(idtFilePtr);

    // create the time history header stored in a character array
    getTimeHistHeader_full(&idt_ctrl, &systemParam, &thistLine);
    // write the time history header to stdout and the time history file
    printf("%s",thistLine.c_str());
    fflush(stdout);
    fprintf(thistFilePtr,"%s",thistLine.c_str());
    fflush(thistFilePtr);

    // create the idt file header stored in a character array
    getIdtHeader(&idt_ctrl, &systemParam, &thistLine);
    // write the idt header to the idt file
    fprintf(idtFilePtr,"%s",thistLine.c_str());
    fflush(idtFilePtr);
  }

  tmax    = idt_ctrl.getStopTime();
  dtprint = idt_ctrl.getPrintTime();

  std::vector<std::string> idt_lines;
  std::vector<std::string> thist_lines;
  int local_idx = 0;

  //Each rank starts at correct offset
  for(j = 0; j < mpi_rank; ++j)
    {
      idt_ctrl.incrementRunId();
    }
  for(j=mpi_rank; j<idt_ctrl.getRunTotal(); j+=mpi_size, local_idx++)
    {
      idt_lines.push_back(std::string(""));
      thist_lines.push_back(std::string(""));

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

      // reset the ignition delay thresholds
      double initTemp = idt_ctrl.getInitTemp();
      std::vector<double> tempDeltas = idt_ctrl.getTemperatureDeltas();
      for(k=0; k < systemParam.numTempRoots; ++k) {
         systemParam.tempRoots[k] = tempDeltas[k]+initTemp;
      }
      idts.clear();
      idts.resize(systemParam.numTempRoots);

      // reset the preconditioner threshold
      change_JsparseThresh(systemParam.sparseMtx,idt_ctrl.getThresh());

      // reset the mass fraction density and relative volume
      idt_ctrl.getInitMassFrac(massFracPtr);
      systemParam.Dens=idt_ctrl.getDensity();
      systemParam.invDens=1.0/systemParam.Dens;
      NV_Ith_S(systemState,nSpc) = initTemp/idt_ctrl.getRefTemp();

      // Save initial state for TLA sensitivity
      for(k=0; k<nState; k++) {
        systemStateCheckpoint[k] = NV_Ith_S(systemState,k);
      }
      timeCheckpoint[0] = 0.0;

      // reset the time
      tcurr=0.0;
      tnext=std::min(dtprint,tmax);

      // reinitialize cvode
      flag = CVodeReInit(cvode_mem, tcurr, systemState);
      isBadStep = 0;

      getTimeHistLine_full(tcurr,NV_DATA_S(systemState),&idt_ctrl,
                           &systemParam,0.0,&thistLine);
#ifdef ZERORK_MPI
      thist_lines[local_idx] += std::string(thistLine);
#else
      printf("%s",thistLine.c_str());
      fflush(stdout);
      fprintf(thistFilePtr,"%s",thistLine.c_str());
      fflush(thistFilePtr);
#endif
      startTime=getHighResolutionTime();

      int count=1;//0;//1;
      int skip=10;
      int count_step=0;
      //double T_old=NV_Ith_S(systemState,nSpc); double t_old=0.0;
      while(tcurr<tmax)
      {
        // Force onestep for sensitivity analysis checkpointing
        flag = CVode(cvode_mem, tmax, systemState, &tcurr, CV_ONE_STEP);
        /*
        if(idt_ctrl.oneStep()) {
          flag = CVode(cvode_mem, tmax, systemState, &tcurr, CV_ONE_STEP);
        }
        else {
          flag = CVode(cvode_mem, tnext, systemState, &tcurr, CV_NORMAL);
        }
        */

        if (check_flag(&flag, "CVODE ERROR", 1)) {

          isBadStep = 1;
          did_cvodeFail++;
          if(did_cvodeFail <= idt_ctrl.getMaxPrimaryCvodeFails()) {

            printf("WARNING: attempting CVodeReInit at t=%.18g [s]\n",tcurr);
            printf("         after cvode failure count = %d for reactor %d\n",
                   did_cvodeFail,j);

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
            printf("         after cvode failure count = %d for reactor %d\n",
                   did_cvodeFail,j);

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
                             &thistLine);

        // Save every "skip" step for TLA sensitivity
        count_step++;
        if(count_step%skip==0) {
          count++;
          count_step=0;
        }
        for(k=0; k<nState; k++) {
          systemStateCheckpoint[count*nState + k] = NV_Ith_S(systemState,k);
        }
        timeCheckpoint[count] = tcurr;
        /**/
        /*
        // Save step when dt or dT exceed some thresholds
        if( (tcurr-t_old > 1.0e-3) or (NV_Ith_S(systemState,nSpc)-T_old>0.002)) {
          count += 1;
          for(k=0; k<nState; k++) {
            systemStateCheckpoint[count*nState + k] = NV_Ith_S(systemState,k);
          }
          timeCheckpoint[count] = tcurr;
          t_old = tcurr;
          T_old = NV_Ith_S(systemState,nSpc);
        }
        */

        // Save dT/dt to evaluate sensitivity coefficient
        dTdt = systemParam.dTemp_dt;

#ifdef ZERORK_MPI
        thist_lines[local_idx] += std::string(thistLine);
#else
        printf("%s",thistLine.c_str());
        fflush(stdout);
        fprintf(thistFilePtr,"%s",thistLine.c_str());
        fflush(thistFilePtr);
#endif

        if(flag == CV_ROOT_RETURN) {
          std::vector<int> roots_found(systemParam.numTempRoots);
          flag = CVodeGetRootInfo(cvode_mem, &roots_found[0]);
          if(flag != CV_SUCCESS) {
            printf("ERROR: CVodeGetRootInfo failed.\n");
            exit(1);
          }
          for(k=0; k < systemParam.numTempRoots; ++k) {
            if(roots_found[k] != 0) { //both rising and falling are kept
              idts[k].push_back(tcurr);
            }
          }
          if(idts[systemParam.numTempRoots-1].size() > 0 &&
             idt_ctrl.continueAfterIDT() == 0) {
            tcurr=2.0*tmax; // advance time past the stopping criteria
          }
        }
        else if(!idt_ctrl.oneStep() && isBadStep==0)
        {tnext+=dtprint;}
      }

      std::vector<double> idts_out(systemParam.numTempRoots,tmax);
      for(k=0; k<systemParam.numTempRoots; ++k) {
         if(idts[k].size() == 0) {
           idts[k].push_back(tmax);
         }
         if(idt_ctrl.getIDTMethod() == 0) {
           idts_out[k] = idts[k][0];
         } else if(idt_ctrl.getIDTMethod() == 1) {
           idts_out[k] = idts[k].back();
         } else { // assume 2
           idts_out[k] = std::accumulate(idts[k].begin(),idts[k].end(),0.0)/idts[k].size();
         }
      }

      // Solve TLA sensitivity ODE
      printf("# Solving TLA sensitivity ODE with %d time checkpoints\n", count-1);
      // TODO: compute sensitivity for multiple temperature deltas?
      if(systemParam.numTempRoots > 1)
        printf("# WARNING: only computing sensitivity analysis for the largest temperature delta\n");
      double tla_start = getHighResolutionTime();
      int nReactions = systemParam.mech->getNumReactions();
      double dt, gamma;
      std::vector<double> Sij, Sij_old, Jz, Z;
      Sij.assign(nReactions*nState,0.0);
      Sij_old.assign(nReactions*nState,0.0);
      Jz.assign(nState,0.0);
      Z.assign(nReactions*nState, 0.0);
      N_Vector fy;
      fy = N_VNew_Serial(nState);
      booleantype jcurPtr;

      // Set sparse threshold to zero to use most exact J
      change_JsparseThresh(systemParam.sparseMtx,0.0);

      // Get y^0, fy^0, and Sij^0
      for(k=0; k<nState; k++) {NV_Ith_S(systemState,k) = systemStateCheckpoint[k];}
      flag = const_vol_wsr_Sij(timeCheckpoint[0], systemState,
                               fy, &Sij[0], &systemParam);

      // Get J^0
      int nsize = systemParam.sparseMtx->nSize;
      N_Vector tmp1, tmp2, tmp3;
      tmp1 = N_VNew_Serial(nsize);
      tmp2 = N_VNew_Serial(nsize);
      tmp3 = N_VNew_Serial(nsize);
      setupJacobianSparse(timeCheckpoint[0], systemState, fy,
                          &systemParam, tmp1, tmp2, tmp3);
      N_VDestroy_Serial(tmp1);
      N_VDestroy_Serial(tmp2);
      N_VDestroy_Serial(tmp3);

      // Time loop
      for(i=0; i<count-1; i++) {
        printf("# Time checkpoint: %d/%d\n",i+1,count-1);
        // Get dt_n+1
        dt = timeCheckpoint[i+1]-timeCheckpoint[i];

        // Get Sij^n
        for(k=0; k<nReactions*nState; k++) {Sij_old[k] = Sij[k];}

        // Get J^n  (copy sparseMtx to sparseMtxOld)
        copy_Jsparse(systemParam.sparseMtx, systemParam.sparseMtxOld);

        // Get y^n+1, fy^n+1, and Sij^n+1
        for(k=0; k<nState; k++) {
          NV_Ith_S(systemState,k) = systemStateCheckpoint[(i+1)*nState + k];
        }
        flag = const_vol_wsr_Sij(timeCheckpoint[i+1], systemState,
                                 fy, &Sij[0], &systemParam);

        // Get J^n+1 and factorize (I-gamma*J^n+1)
        gamma = 0.5*dt;
        flag = jac_full_prec_setup(timeCheckpoint[i+1], systemState, fy,
                                   false, &jcurPtr, gamma, &systemParam);

        // Reaction loop
        for(k=0; k<nReactions; k++) {

          // Get J^n*z^n
          flag = sparse_oldjac_v(&Z[k*nState], &Jz[0], &systemParam);

          // Compute RHS (Z^n + 1/2*dt_n+1*(J^n*Z^n + S^n+1 + S^n))
          for(int ii=0; ii<nState; ii++) {
            Z[k*nState+ii] += 0.5*dt*(Jz[ii] + Sij[k*nState+ii] +
                                      Sij_old[k*nState+ii]);
          }

          // Solve linear system (I-gamma*J)*Z^n+1 = RHS to get Z^n+1
          backSolve(&Z[k*nState], &systemParam);

        } // for k<nReactions
      } // for i<count-1

      N_VDestroy_Serial(fy);

      // Get IDT sensitivity coefficients
      rxnSens_t *rxnSensList;
      rxnSensList = new rxnSens_t[nReactions];
      for(k=0; k<nReactions; k++) {
        // S_k = -1/tau * (dT/dxi)/(dT/dt)
        rxnSensList[k].relSens = -1/timeCheckpoint[count]*Z[k*nState + nSpc]/dTdt;
        rxnSensList[k].rxnId = k;
      }

      // Sort sensitivity coefficients in descending order
      qsort((void *)&rxnSensList[0], nReactions, sizeof(rxnSens_t), compare_rxnSens_t);

      // print to sens coeffs to file -- separate file for each run
      char sensFileName[32];
      snprintf(sensFileName, 32, "sensAnalysis_IDT_run%03d.txt",idt_ctrl.getRunId());
      FILE * sensFile;
      sensFile = fopen(sensFileName,"w");
      fprintf(sensFile,"#reac_id  Srel       Reaction\n");
      for(int k=0; k<nReactions; k++) {
        fprintf(sensFile,"%8d  %16.7e  %s\n",
                rxnSensList[k].rxnId+1,//to start at 1 instead of 0
                rxnSensList[k].relSens,
                systemParam.mech->getReactionName(rxnSensList[k].rxnId));
      }
      fclose(sensFile);

      double TLA_time = getHighResolutionTime() - tla_start;
      printf("# Sensitivity analysis took %5.3e seconds\n", TLA_time);
      // Done solving TLA sensitivity

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
#ifdef ZERORK_MPI
      thist_lines[local_idx] += std::string("\n");
#else
      printf("\n");
      fprintf(thistFilePtr,"\n");
#endif

      idt_ctrl.getInitMoleFrac(&moleFracInit[0]);
      getIdtLine(&moleFracInit[0],&idt_ctrl,idts_out,simTime,did_cvodeFail, &thistLine);
#ifdef ZERORK_MPI
      idt_lines[local_idx] = std::string(thistLine);
#else
      fprintf(idtFilePtr,"%s",thistLine.c_str());
      fflush(idtFilePtr);
#endif

      avgTime[idt_ctrl.getThreshId()]+=simTime;
      if(simTime > maxTime[idt_ctrl.getThreshId()])
	{maxTime[idt_ctrl.getThreshId()]=simTime;}
      if(simTime < minTime[idt_ctrl.getThreshId()])
	{minTime[idt_ctrl.getThreshId()]=simTime;}

      for(int jj = 0; jj < mpi_size; ++jj) {
        idt_ctrl.incrementRunId();
      }
    }

#ifdef ZERORK_MPI
  for(j = 0; j < (idt_ctrl.getRunTotal()+mpi_size-1)/mpi_size; ++j)
  {
    for(k = 0; k < mpi_size; ++k)
    {
      int current_run = mpi_size*j+k;
      if(current_run >= idt_ctrl.getRunTotal())
      {
        break;
      }
      std::string current_thist_lines;
      std::string current_idt_lines;
      if(k != 0) {
        mpi_sendrecv_string(k, 0, thist_lines[j], &current_thist_lines, MPI_COMM_WORLD);
        mpi_sendrecv_string(k, 0, idt_lines[j], &current_idt_lines, MPI_COMM_WORLD);
      } else {
        if(mpi_rank == 0) {
          current_thist_lines = thist_lines[j];
          current_idt_lines = idt_lines[j];
        }
      }
      if(mpi_rank == 0)
      {
          printf("%s",current_thist_lines.c_str());
          fflush(stdout);
          fprintf(thistFilePtr,"%s",current_thist_lines.c_str());
          fflush(thistFilePtr);

          fprintf(idtFilePtr,"%s",current_idt_lines.c_str());
          fflush(idtFilePtr);
      }
    }
  }

  MPI_Reduce(&maxTime[0],&global_maxTime[0],maxTime.size(),MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&minTime[0],&global_minTime[0],minTime.size(),MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&avgTime[0],&global_avgTime[0],avgTime.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  maxTime = global_maxTime; //over-writing to save mpi logic during printing
  minTime = global_minTime;
  avgTime = global_avgTime;
#endif

  if(mpi_rank == 0) {
    // report the threshold statistics
    fprintf(idtFilePtr,"\n\n");
    fprintf(idtFilePtr,"# -----------------------------------------------------------------------------\n");
    fprintf(idtFilePtr,"# Column 1: [-] preconditioner threshold\n");
    fprintf(idtFilePtr,"# Column 2: [s] average simulation time\n");
    fprintf(idtFilePtr,"# Column 3: [s] minimum simulation time\n");
    fprintf(idtFilePtr,"# Column 4: [s] maximum simulation time\n");
    nRunsPerThresh=(double)(idt_ctrl.getRunTotal()/
          		  idt_ctrl.getNumThreshRuns());
    for(j=0; j<idt_ctrl.getNumThreshRuns(); j++) {
      avgTime[j]/=nRunsPerThresh;
      fprintf(idtFilePtr,"%14.7e  %14.7e  %14.7e  %14.7e\n",
              idt_ctrl.getThreshFromId(j),avgTime[j],minTime[j],maxTime[j]);
    }
    fflush(idtFilePtr);
    fclose(idtFilePtr);
    fclose(thistFilePtr);
  }

  // clean up procedures
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
#ifdef ZERORK_MPI
  MPI_Init(&argc, &argv);
#endif
  cvReactor(argc, argv);
#ifdef ZERORK_MPI
  MPI_Finalize();
#endif
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
                          std::string *thistLine)
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
  if(nsteps == 0) {
    nniters = 0; //Fix bug in CVode that doesn't reset this on ReInit
  }
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

  std::stringstream ss;

  ss << std::setw(8);
  ss << idt_ctrl->getRunId();
  int width = 16;
  ss << std::setprecision(7);
  ss << std::scientific;
  if(idt_ctrl->longOutput()) {
    width = 26;
    ss << std::setprecision(16);
  }
  ss << std::setw(width) << currTime
     << std::setw(width) << currTemp
     << std::setw(width) << currPres
     << std::setw(width) << 1.0/cvp->invDens
     << std::setw(width) << molWt
     << std::setw(width) << cvMixPerMass
     << std::setw(width) << intEnergyMixPerMass-intEnergyMixPerMassInit;
  for(j=0; j<nTrack; j++) {
      ss << std::setw(width) << currMoleFrac[idt_ctrl->getTrackSpeciesId(j)];
   }
  if(idt_ctrl->printNetProductionRates()) {
    for(j=0; j<nTrack; j++) {
      ss << std::setw(width) << cvp->netProd[idt_ctrl->getTrackSpeciesId(j)];
    }
  }
  if(idt_ctrl->printNetRatesOfProgress()) {
    int nRxn = cvp->mech->getNumReactions();
    int nStep=idt_ctrl->getNumSteps();
    int fwdStepIdx, revStepIdx;
    double netROP;
    for(j=0; j<nRxn; j++) {
      cvp->mech->getStepIdxOfRxn(j, &fwdStepIdx, &revStepIdx);
      netROP = cvp->fwdROP[fwdStepIdx];
      if(revStepIdx > 0) {
         netROP -= cvp->fwdROP[revStepIdx];
      }
      ss << std::setw(width) << netROP;
    }
  }

  ss << std::setw(8) << nsteps
     << std::setw(8) << netfails
     << std::setw(8) << nniters
     << std::setw(8) << nncfails
     << std::setw(8) << cvp->nJacSetup
     << std::setw(8) << cvp->nColPerm
     << std::setw(8) << cvp->nJacFactor
     << std::setw(8) << cvp->nBackSolve
     << std::setw(8) << cvp->nFunc
     << std::setw(8) << cvp->sparseMtx->nNonZero
     << std::setw(8) << cvp->sparseMtx->reduceNNZ
     << std::setw(8) << cvp->sparseMtx->LUnnz;

  ss << std::setprecision(7);
  ss << std::setw(16) << cvp->jacSetupTime
     << std::setw(16) << cvp->colPermTime
     << std::setw(16) << cvp->jacFactorTime
     << std::setw(16) << cvp->backsolveTime
     << std::setw(16) << cvp->funcTime
     << std::setw(16) << wallTime
     << std::endl;

  *thistLine = ss.str();
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
                            std::string *thistHead)
{
  int j;
  int nSpc   = cvp->mech->getNumSpecies();
  int nTrack = idt_ctrl->getNumTrackSpecies();
  int colNum;
  std::stringstream ss;

  ss << "# Column  1: [#]       run id\n"
     << "# Column  2: [s]       simulated ODE time\n"
     << "# Column  3: [K]       temperature\n"
     << "# Column  4: [Pa]      pressure\n"
     << "# Column  5: [kg/m^3]  density\n"
     << "# Column  6: [kg/kmol] mixture molecular weight\n"
     << "# Column  7: [J/kg/K]  mixture specific heat (cosnt vol)\n"
     << "# Column  8: [J/kg]    energy released\n";

  colNum=9;
  ss << std::setw(2);
  for(j=0; j<nTrack; j++) {
    ss <<"# Column " << colNum++ << ": [-]      mole fraction of "
       << cvp->mech->getSpeciesName(idt_ctrl->getTrackSpeciesId(j))
       << std::endl;
  }
  if(idt_ctrl->printNetProductionRates()) {
    for(j=0; j<nTrack; j++) {
      ss <<"# Column " << colNum++ << ": [kmol/m^3/s]      net production rate of "
         << cvp->mech->getSpeciesName(idt_ctrl->getTrackSpeciesId(j))
         << std::endl;
    }
  }
  if(idt_ctrl->printNetRatesOfProgress()) {
    int nRxn = cvp->mech->getNumReactions();
    for(j=0; j<nRxn; j++) {
      ss <<"# Column " << colNum++ << ": [kmol/m^3/s]      net rate of progress for reaction "
         << j << std::endl;
    }
  }
  ss << "# Column " << colNum++ << ": [#]      CVodeGetNumSteps(...)\n"
     << "# Column " << colNum++ << ": [#]      CVodeGetNumErrTestFails(...)\n"
     << "# Column " << colNum++ << ": [#]      CVodeGetNumNonlinSolvIters(...)\n"
     << "# Column " << colNum++ << ": [#]      CVodeGetNumNonlinSolvConvFails(...)\n"
     << "# Column " << colNum++ << ": [#]      calls to the precond. setup\n"
     << "# Column " << colNum++ << ": [#]      calls to the column permutation\n"
     << "# Column " << colNum++ << ": [#]      calls to the precond. factorization\n"
     << "# Column " << colNum++ << ": [#]      calls to the precond. back solve\n"
     << "# Column " << colNum++ << ": [#]      calls to the RHS function\n"
     << "# Column " << colNum++ << ": [#]      max non-zero terms in precond.\n"
     << "# Column " << colNum++ << ": [#]      current non-zero terms in precond.\n"
     << "# Column " << colNum++ << ": [#]      current non-zero terms in LU factors\n"
     << "# Column " << colNum++ << ": [s]      wall clock time - precond. setup\n"
     << "# Column " << colNum++ << ": [s]      wall clock time - column permut.\n"
     << "# Column " << colNum++ << ": [s]      wall clock time - precond. factor.\n"
     << "# Column " << colNum++ << ": [s]      wall clock time - precond. back solve\n"
     << "# Column " << colNum++ << ": [s]      wall clock time - RHS evaluation\n"
     << "# Column " << colNum++ << ": [s]      wall clock time - total simulation\n"
     << "# -----------------------------------------------------------------------------\n"
     << "# run id";

  int width = 16;
  if(idt_ctrl->longOutput()) {
    width = 26;
  }
  ss << std::setw(width) << "sim time"
     << std::setw(width) << "temp"
     << std::setw(width) << "press"
     << std::setw(width) << "dens"
     << std::setw(width) << "mol wt"
     << std::setw(width) << "cv"
     << std::setw(width) << "E release";
  for(j=0; j<nTrack; j++) {
    ss << " mlfrc" << std::setw(width-6) << cvp->mech->getSpeciesName(idt_ctrl->getTrackSpeciesId(j));
  }
  if(idt_ctrl->printNetProductionRates()) {
    for(j=0; j<nTrack; j++) {
      ss << " npr  " << std::setw(width-6) << cvp->mech->getSpeciesName(idt_ctrl->getTrackSpeciesId(j));
    }
  }
  if(idt_ctrl->printNetRatesOfProgress()) {
    int nRxn = cvp->mech->getNumReactions();
    for(j=0; j<nRxn; j++) {
      ss << " rop  " << std::setw(width-6) << j;
    }
  }
  ss << "   nstep   nertf  nnlsit  nnlcvf prsetup ncolprm nfactor backslv  intfcl     nnz    rnnz   lunnz    jac setup tm     col perm tm      jac fac tm    backsolve tm     rhs func tm   simulation tm\n";
 *thistHead = ss.str();
}

// void getSimHeader_full(...)
//
//
void getSimHeader_full(const int inp_argc,
                       char **inp_argv,
                       const idt_sweep_params *idt_ctrl,
                       const cv_param *cvp,
                       std::string *simHead)
{
  int j,spcId;
  int nSpc   = cvp->mech->getNumSpecies();
  int nSpcTermsJac = cvp->sparseMtx->nNonZero-2*nSpc+1;
  int nTrack = idt_ctrl->getNumTrackSpecies();
  std::stringstream ss;

  // 1. construct the command line
  ss << std::string("# Command line: ") << inp_argv[0];

  for(j=1; j<inp_argc; j++) {
    ss << " " << inp_argv[j];
  }
  ss << "\n# ---------------------------------------------------------------------------\n";

  {
    zerork::ZeroRKConstants zrk_constants;
    // 1.5 construct git commit information
    ss << "# Zero-RK git commit information:" <<
          "\n#   commit id                 [hash] : " <<
          zrk_constants.GetCommitId() <<
          "\n#   commit time-stamp         [date] : " <<
          zrk_constants.GetCommitTimestamp() <<
          "\n#   branch name               [word] : " <<
          zrk_constants.GetBranchName() <<
          "\n#   exponential type          [word] : " <<
          zrk_constants.GetExpType() <<
          "\n# ---------------------------------------------------------------------------\n";
  }


  // 2. construct the mechanism information
  ss << "# Mechanism constants:" <<
        "\n#   mechanism definition      [path] : " <<
        idt_ctrl->getMechFileName() <<
        "\n#   thermodynamics definition [path] : " <<
        idt_ctrl->getThermFileName() <<
        "\n#   parser log definition     [path] : " <<
        idt_ctrl->getLogFileName() <<
        "\n#   number of species            [#] : " <<
        std::to_string(nSpc) <<
        "\n#   number of reactions          [#] : " <<
        std::to_string(cvp->mech->getNumReactions()) <<
        "\n#   number of 1-way steps        [#] : " <<
        std::to_string(cvp->mech->getNumSteps()) <<
        "\n#   max non-zeros in species J   [#] : " <<
        std::to_string(nSpcTermsJac) <<
        "\n#   max non-zeros in const vol J [#] : " <<
        std::to_string(cvp->sparseMtx->nNonZero) <<
        "\n# ---------------------------------------------------------------------------\n";

  // 3. construct the composition information
  ss << std::setprecision(4) << std::fixed;
  ss << "# Fuel composition (moles of species)/(moles of fuel):" << std::endl;
  for(j=0; j<idt_ctrl->getNumFuelSpecies(); j++) {
    spcId = idt_ctrl->getFuelSpeciesId(j);
    ss << "#  " << std::setw(16) << cvp->mech->getSpeciesName(spcId);
    ss << ": " << idt_ctrl->getFuelMoleFrac(spcId) << std::endl;
  }
  ss << "# Oxid(izer) composition (moles of species)/(moles of oxid):" << std::endl;
  for(j=0; j<idt_ctrl->getNumOxidSpecies(); j++) {
    spcId = idt_ctrl->getOxidSpeciesId(j);
    ss << "#  " << std::setw(16) << cvp->mech->getSpeciesName(spcId)
       << ": " << idt_ctrl->getOxidMoleFrac(spcId) << std::endl;
  }
  ss << std::setw(0);

  ss << "# Global stoichiometry constants:" << std::endl
     << "#   fuel (excess mols of atomic-O)/(mols of fuel) : "
     << idt_ctrl->getFuelOxyBalance() << std::endl
     << "#   oxid (excess mols of atomic-O)/(mols of oxid) : "
     << idt_ctrl->getOxidOxyBalance() << std::endl
     << "#   (mols of oxid)/(mols of fuel) at phi=1        : "
     << idt_ctrl->getMoleOxidStoicPerFuel() << std::endl
     << "# ---------------------------------------------------------------------------" << std::endl;

  // 4. construct the simulation information
  ss << "# Simulation constants:" << std::endl
     << std::setprecision(4) << std::scientific
     << "#   ode print time               [s] : "
     << idt_ctrl->getPrintTime() << std::endl
     << "#   ode max time                 [s] : "
     << idt_ctrl->getStopTime() << std::endl
     << "#   continue after IDT found   [y/n] : "
     << ((idt_ctrl->continueAfterIDT() != 0) ? "yes" : "no") << std::endl
     << std::setprecision(2) << std::fixed
     << "#   reference temperature        [K] : "
     << idt_ctrl->getRefTemp() << std::endl;
  for(j=0; j<idt_ctrl->getTemperatureDeltas().size(); j++) {
      ss << "#   ignition delay delta Temp[" << j << "]    [K] : "
         << idt_ctrl->getTemperatureDeltas()[j] << std::endl;
  }
  ss << "#   total number of runs         [#] : "
     << idt_ctrl->getRunTotal() << std::endl
     << "# ---------------------------------------------------------------------------" << std::endl;

  // 5. construct the solver information
  ss << std::setprecision(4) << std::scientific;
  ss << "# Solver constants:" << std::endl
     << "#   ode relative tolerance            [-] : "
     << idt_ctrl->getRTol() << std::endl
     << "#   ode absolute tolerance            [-] : "
     << idt_ctrl->getATol() << std::endl
     << "#   ode max internal timestep         [s] : "
     << idt_ctrl->getMaxInternalDt() << std::endl
     << "#   ode max internal steps            [#] : "
     << idt_ctrl->getMaxInternalSteps() << std::endl
     << "#   linear solver rel error           [-] : "
     << idt_ctrl->getEpsLin() << std::endl
     << "#   non-linear solver rel error       [-] : "
     << idt_ctrl->getNlConvCoeff() << std::endl
     << "#   cvode primary failures allowed    [#] : "
     << idt_ctrl->getMaxPrimaryCvodeFails() << std::endl
     << "#   cvode secondary failures allowed  [#] : "
     << idt_ctrl->getMaxSecondaryCvodeFails() << std::endl
     << "#   secondary failure P threshold     [-] : "
     << idt_ctrl->getSafetyThreshold() << std::endl
     << "#   max Krylov dimension              [#] : "
     << idt_ctrl->getKrylovDim() << std::endl
     << "#   use fake P update               [y/n] : "
     << (idt_ctrl->getUpdate() ? "yes" : "no") << std::endl
     << "#   use incomplete LU factorization [y/n] : "
     << (idt_ctrl->getILU() ? "yes" : "no") << std::endl
     << "#   partial pivot threshold           [-] : "
     << idt_ctrl->getPartialPivotThresh() << std::endl
     << "#   P threshold type code             [#] : "
     << precondThresholdName[idt_ctrl->getThreshType()-1] << std::endl
     << "#   permutation type code             [#] : "
     << permutationTypeName[idt_ctrl->getPermutationType()-1] << std::endl
     << "# ---------------------------------------------------------------------------" << std::endl;

  *simHead = ss.str();
}

void getIdtLine(const double initMoleFrac[],
               const idt_sweep_params *idt_ctrl,
               const std::vector<double> idts,
               const double runTime,
               const int nCvodeFails,
               std::string *idtLine)
{
  int j;
  int nTrack = idt_ctrl->getNumTrackSpecies();
  std::stringstream ss;
  ss << std::setprecision(16);
  ss << std::scientific;

  ss << std::setw(8) << idt_ctrl->getRunId()
     << std::setw(25) << idt_ctrl->getInitTemp()
     << std::setw(25) << idt_ctrl->getInitPres()
     << std::setw(25) << idt_ctrl->getPhi()
     << std::setw(25) << idt_ctrl->getEgr()
     << std::setw(25) << idt_ctrl->getThresh()
     << std::setw(25) << idt_ctrl->getDensity();
  for(j=0; j<idts.size(); ++j) {
      ss << std::setw(25) << idts[j];
  }
  ss << std::setw(25) << runTime
     << std::setw(16) << nCvodeFails;

  for(j=0; j<nTrack; j++) {
    ss << std::setw(25) << initMoleFrac[idt_ctrl->getTrackSpeciesId(j)];
  }
  ss << std::endl;
  *idtLine = ss.str();
}

void getIdtHeader(const idt_sweep_params *idt_ctrl,
                 const cv_param *cvp,
                 std::string *idtHeader)
{
  int j;
  int nTrack = idt_ctrl->getNumTrackSpecies();
  int numTempDeltas = idt_ctrl->getTemperatureDeltas().size();
  std::stringstream ss;

  ss << "# Column  1: [#]      run id" << std::endl
     << "# Column  2: [K]      initial temperature" << std::endl
     << "# Column  3: [Pa]     initial pressure" << std::endl
     << "# Column  4: [-]      equivalence ratio" << std::endl
     << "# Column  5: [-]      EGR fraction" << std::endl
     << "# Column  6: [-]      preconditioner threshold" << std::endl
     << "# Column  7: [kg/m^3] initial density" << std::endl;
  for(j=0; j<numTempDeltas; ++j) {
     ss << "# Column " << std::setw(2) << j+8
        <<": [s]      ignition delay time (delta T = "
        << idt_ctrl->getTemperatureDeltas()[j]
        <<" K)" << std::endl;
  }
  ss << "# Column " << std::setw(2) << 8+numTempDeltas
     << ": [s]      execution wall clock time" << std::endl;
  ss << "# Column " << std::setw(2) << 9+numTempDeltas
     << ": [#]      number of cvode errors" << std::endl;

  for(j=0; j<nTrack; j++) {
    ss << "# Column " << j+numTempDeltas+10 << ": [-]      initial mole frac of "
       << cvp->mech->getSpeciesName(idt_ctrl->getTrackSpeciesId(j)) << std::endl;
  }
  ss << "# ----------------------------------------------------------------------------\n";

  ss << std::setw(7) << "# run id |"
     << std::setw(25) << "Temp [K] |"
     << std::setw(25) << "Pres [Pa] |"
     << std::setw(25) << "phi [-] |"
     << std::setw(25) << "EGR frac [-] |"
     << std::setw(25) << "P thresh [-] |"
     << std::setw(25) << "dens [kg/m3] |";
  for(j=0; j<numTempDeltas; ++j) {
      ss << std::setw(13) << idt_ctrl->getTemperatureDeltas()[j]
         << " K IDT [s] |";
  }
  ss << std::setw(25) << "sim time [s] |"
     << std::setw(16) << "cvode err [#] |";

  for(j=0; j<nTrack; j++) {
    ss << " mlfrc " << std::setw(16)
       << cvp->mech->getSpeciesName(idt_ctrl->getTrackSpeciesId(j))
       << " |";
  }
  ss << std::endl;
  *idtHeader = ss.str();
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


#ifdef ZERORK_MPI
static void mpi_sendrecv_string(int from, int to, std::string& in, std::string* out, MPI_Comm comm)
{
  int mpi_rank = 0;
  int mpi_size = 1;
  MPI_Comm_rank(comm,&mpi_rank);
  MPI_Comm_size(comm,&mpi_size);
  assert(mpi_size > 1);
  assert(to < mpi_size);
  assert(from < mpi_size);

  if(mpi_rank == from) {
    MPI_Send(const_cast<void*>((void*)in.c_str()), in.size()+1, MPI_CHAR, to, 0, comm);
  }
  if(mpi_rank == to) {
    MPI_Status status;
    MPI_Probe(from, 0, comm, &status);
    int len = 0;
    MPI_Get_count(&status, MPI_CHAR, &len);
    char *buf = new char[len];
    MPI_Recv(buf, len, MPI_CHAR, from, 0, comm, &status);
    *out = buf;
    delete [] buf;
  }
  MPI_Barrier(comm);
}
#endif


int compare_rxnSens_t(const void *A, const void *B)
{
  rxnSens_t *Aptr =(rxnSens_t *)A;
  rxnSens_t *Bptr =(rxnSens_t *)B;
  if(fabs(Aptr->relSens) < fabs(Bptr->relSens)) {
    return 1;
  }
  else if (fabs(Aptr->relSens) > fabs(Bptr->relSens)) {
    return -1;
  }
  return 0;
}
