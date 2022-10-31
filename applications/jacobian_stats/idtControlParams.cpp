#include <string.h>

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

#include "utility_funcs.h"
#include "ode_funcs.h"
#include "matrix_funcs.h"
#include "idtControlParams.h"
#include "mechanism_stats.h"


idtControlParams::idtControlParams(char *inpFile, int printLevel)
{
  BasicReactorIFP parser(inpFile);

  // -----------------------------------------------------------------------
  // zerork mechanism setup:
  mechFile    = parser.mechFile();
  thermFile   = parser.thermFile();
  mechLogFile = parser.mechLogFile();
  mechStatFile = parser.mechStatFile();
  jacobianStatFile = parser.jacobianStatFile();
  jacobianRawFile  = parser.jacobianRawFile();
  integratorProbFile = parser.integratorProbFile();
  outFile     = parser.outFile();
  printDeltaTime = parser.printDeltaTime();

  if(printLevel == 0) {
    mech = new zerork::mechanism(mechFile.c_str(),
                              thermFile.c_str(),
                              "");
  } else {
    mech = new zerork::mechanism(mechFile.c_str(),
                              thermFile.c_str(),
                              mechLogFile.c_str());
  }
  // set local copies of basic mechanism size constants.
  nSpc  = mech->getNumSpecies();
  nRxn  = mech->getNumReactions();
  nStep = mech->getNumSteps();

  // set the temperature scan parameters for mechanism checks
  num_scans = parser.num_scans();
  scan_temperature_ref = parser.scan_temperature_ref();
  min_scan_temperature = parser.min_scan_temperature();
  max_scan_temperature = parser.max_scan_temperature();


  // get the AFactor perturbation settings.
  //doBothDir = parser.doBothDir();
  //AFactorMultiplier = parser.AFactorMultiplier();
  //if(AFactorMultiplier == 1.0) {
  //  printf("ERROR: AFactorMultiper can not equal one.\n");
  //  exit(-1);
  //}

  ropMultiplier = new double[nStep];
  clearAllROPMultiplier();  // set ropMultiplier[:] = 1.0


  setupZeroRKParams(parser);
  setupCVode(parser);       // includes allocating N_Vector systemState
  setupIdtState(parser);
  setupCVodeUserParams(parser);
  resetIdtSim(); // sets up the system state

}

idtControlParams::~idtControlParams(void)
{
  freeCVodeUserParams();

  delete [] ropMultiplier;
  delete [] fuelMoleFrac;
  delete [] oxidMoleFrac;
  delete [] initMoleFrac;
  delete [] initMassFrac;

  delete mech;

  N_VDestroy_Serial(systemState);
  // Free integrator memory
  CVodeFree(&cvodeCtrl.cvodeMemPtr);
#if defined SUNDIALS3 || defined SUNDIALS4
  SUNLinSolFree(cvodeCtrl.LS);
#endif
#if defined SUNDIALS4
  SUNNonlinSolFree(cvodeCtrl.NLS);
#endif
}

// Reset the sundials' NVector for systemState, re-initialized cvode
// to time zero, and set all the solver counters to zero
void idtControlParams::resetIdtSim(void)
{
  int j;
  int flag;

  for(j=0; j<nSpc; j++) {
    NV_Ith_S(systemState,j)=initMassFrac[j];
  }
  NV_Ith_S(systemState,nSpc)=initTemp/refTemp;
  flag = CVodeReInit(cvodeCtrl.cvodeMemPtr, 0.0, systemState);
  if (check_flag(&flag, "CVodeReInit", 1)) {
    exit(-1);
  }
  odeUserParams.sparseMtx->reduceNNZ=0; // won't be set until the first J
  odeUserParams.sparseMtx->LUnnz = 0; // won't be set until the first J
  odeUserParams.sparseMtx->fillFactor= 100.; // won't be set until the first J
  odeUserParams.prevNumErrTestFails = 0;
  odeUserParams.nFunc=0;
  odeUserParams.nJacSetup=0;
  odeUserParams.nJacFactor=0;
  odeUserParams.nBackSolve=0;
  odeUserParams.nJacRescale=0;
  odeUserParams.nColPerm=0;
  odeUserParams.colPermTime   = 0.0;
  odeUserParams.jacFactorTime = 0.0;
  odeUserParams.backsolveTime = 0.0;
  odeUserParams.jacSetupTime = 0.0;
  odeUserParams.funcTime = 0.0;

  num_cvode_warnings = 0;
  num_cvode_errors   = 0;
}

void idtControlParams::clearAllROPMultiplier(void)
{
  for(int j=0; j<nStep; j++) {
    ropMultiplier[j]=1.0;
  }
}

void idtControlParams::setROPMultiplierOfRxn(const int rxnId,
                                              const bool doDivision)
{
  double Amult=AFactorMultiplier;
  int fwdId = mech->getStepIdxOfRxn(rxnId,1);
  int revId = mech->getStepIdxOfRxn(rxnId,-1);

  if(doDivision) {
    Amult = 1.0/Amult;
  }
  ropMultiplier[fwdId]*=Amult;
  if(revId >= 0 && revId < nStep) {
    ropMultiplier[revId]*=Amult;
  }
}
void idtControlParams::unsetROPMultiplierOfRxn(const int rxnId)
{
  int fwdId = mech->getStepIdxOfRxn(rxnId,1);
  int revId = mech->getStepIdxOfRxn(rxnId,-1);

  ropMultiplier[fwdId]=1.0;
  if(revId >= 0 && revId < nStep) {
    ropMultiplier[revId]=1.0;
  }
}

void idtControlParams::setupCVode(BasicReactorIFP &parser)
{
  int flag;

  // set cvode control parameters
  cvodeCtrl.nDims             = nSpc+1;
  cvodeCtrl.nRoots            = parser.idtTemps().size();
  cvodeCtrl.maxSteps          = parser.maxSteps();
  cvodeCtrl.krylovDim         = parser.krylovDim();
  cvodeCtrl.maxNumNonLinIters = parser.maxNumNonLinIters();
  cvodeCtrl.relTol            = parser.relTol();
  cvodeCtrl.absTol            = parser.absTol();
  cvodeCtrl.maxDtInternal     = parser.maxDtInternal();
  cvodeCtrl.cvEpsLin          = parser.cvEpsLin();
  cvodeCtrl.cvNlConvCoeff     = parser.cvNlConvCoeff();
  cvodeCtrl.maxTime           = parser.maxTime();

  // NVector for cvode system state
  systemState = N_VNew_Serial(cvodeCtrl.nDims);

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration */
#ifdef SUNDIALS4
  cvodeCtrl.cvodeMemPtr = CVodeCreate(CV_BDF);
#else
  cvodeCtrl.cvodeMemPtr = CVodeCreate(CV_BDF, CV_NEWTON);
#endif
  if(check_flag((void *)cvodeCtrl.cvodeMemPtr,"CVodeCreate",0)) {
    exit(-1);
  }

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag=CVodeInit(cvodeCtrl.cvodeMemPtr,
                 const_vol_wsr_limiter,
                 0.0,
                 systemState);
  if (check_flag(&flag, "CVodeInit", 1)) {
    exit(-1);
  }

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  flag = CVodeSStolerances(cvodeCtrl.cvodeMemPtr,
			   cvodeCtrl.relTol,
                           cvodeCtrl.absTol);
  if (check_flag(&flag, "CVodeSStolerances", 1)) {
    exit(-1);
  }

  /* Call CVodeRootInit to specify the root function with 1 component */
  flag = CVodeRootInit(cvodeCtrl.cvodeMemPtr,
                       cvodeCtrl.nRoots,
                       tempRootFunc);
  if (check_flag(&flag, "CVodeRootInit", 1)) {
    exit(-1);
  }

  /* Set the pointer to user-defined data */
  flag = CVodeSetUserData(cvodeCtrl.cvodeMemPtr,
                          &odeUserParams);
  if(check_flag(&flag, "CVodeSetUserData", 1)) {
    exit(-1);
  }

  /* Call CVSpgmr to specify the linear solver CVSPGMR
     with left preconditioning and the maximum Krylov dimension maxl */
#ifdef SUNDIALS2
  flag = CVSpgmr(cvodeCtrl.cvodeMemPtr,
                 PREC_LEFT,
                 cvodeCtrl.krylovDim);
  if(check_flag(&flag, "CVSpgmr", 1)) {
    exit(-1);
  }
  /* Call CVSpilsSetGSType to specify the QR factorization */
  flag = CVSpilsSetGSType(cvodeCtrl.cvodeMemPtr,
                          MODIFIED_GS);
  if(check_flag(&flag, "CVSpilsSetGSType", 1)) {
    exit(-1);
  }

  /* Set preconditioner setup and solve routines Precond and PSolve,
     and the pointer to the user-defined block data */
  flag = CVSpilsSetPreconditioner(cvodeCtrl.cvodeMemPtr,
                                  jac_full_prec_setup,
                                  jac_full_prec_solveV3);
  if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) {
    exit(-1);
  }
#elif SUNDIALS3
  cvodeCtrl.LS = SUNSPGMR(systemState, PREC_LEFT, cvodeCtrl.krylovDim);
  flag = CVSpilsSetLinearSolver(cvodeCtrl.cvodeMemPtr, cvodeCtrl.LS);
  if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) exit(-1);

  flag = SUNSPGMRSetGSType(cvodeCtrl.LS, MODIFIED_GS);
  if(check_flag(&flag, "SUNSPGMRSetGSType", 1)) exit(-1);
    /* Set preconditioner setup and solve routines Precond and PSolve,
       and the pointer to the user-defined block data */
  flag = CVSpilsSetPreconditioner(cvodeCtrl.cvodeMemPtr,
                                  jac_full_prec_setup,
                                  jac_full_prec_solveV3);
  if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) exit(-1);
#elif SUNDIALS4
  cvodeCtrl.NLS = SUNNonlinSol_Newton(systemState);
  flag = CVodeSetNonlinearSolver(cvodeCtrl.cvodeMemPtr, cvodeCtrl.NLS);
  if(check_flag(&flag, "CVodeSetNonlinearSolver", 1)) exit(-1);
  cvodeCtrl.LS = SUNLinSol_SPGMR(systemState, PREC_LEFT, cvodeCtrl.krylovDim);
  flag = CVodeSetLinearSolver(cvodeCtrl.cvodeMemPtr, cvodeCtrl.LS, NULL);
  if(check_flag(&flag, "CVodeSetLinearSolver", 1)) exit(-1);

  /* Set preconditioner setup and solve routines Precond and PSolve,
     and the pointer to the user-defined block data */
  flag = CVodeSetPreconditioner(cvodeCtrl.cvodeMemPtr,
                                jac_full_prec_setup,
                                jac_full_prec_solveV3);
  if(check_flag(&flag, "CVodeSetPreconditioner", 1)) exit(-1);
#endif

  //  flag = CVSpilsSetJacTimesVecFn(cvode_mem, sparse_jac_v);
  //  if(check_flag(&flag, "CVSpilsSetJacTimesVecFn", 1)) exit(-1);

  /* Set the maximum number of internal steps per CVode call and the maximum
   * allowable internal steps. */
  flag = CVodeSetMaxNumSteps(cvodeCtrl.cvodeMemPtr,
                             cvodeCtrl.maxSteps);
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) {
    exit(-1);
  }
  flag = CVodeSetMaxStep(cvodeCtrl.cvodeMemPtr,
                         cvodeCtrl.maxDtInternal);
  if (check_flag(&flag, "CVodeSetMaxStep", 1)) {
    exit(-1);
  }

  /* Set the safety factors on the iterative solver tolerances for the
     nonlinear and linear solvers */

  flag = CVodeSetNonlinConvCoef(cvodeCtrl.cvodeMemPtr,
                                cvodeCtrl.cvNlConvCoeff);
  if (check_flag(&flag, "CVodeSetNonlinConvCoef", 1)) {
    exit(-1);
  }

#if defined SUNDIALS2 || defined SUNDIALS3
  flag = CVSpilsSetEpsLin(cvodeCtrl.cvodeMemPtr,
                          cvodeCtrl.cvEpsLin);
#else
  flag = CVodeSetEpsLin(cvodeCtrl.cvodeMemPtr,
                        cvodeCtrl.cvEpsLin);
#endif

}

void idtControlParams::setupCVodeUserParams(BasicReactorIFP &parser)
{
  int j;

  // set ode user parameter constants
  odeUserParams.mechPtr=mech;
  odeUserParams.nSpc = nSpc;
  odeUserParams.Tref = refTemp;
  odeUserParams.Dens = initDens;
  odeUserParams.invDens = 1.0/initDens;
  odeUserParams.minMassFrac = 1.0e-30;
  odeUserParams.sqrtUnitRnd = sqrt(UNIT_ROUNDOFF);

  // allocate ode user parameter memory
  odeUserParams.netProd     = new double[nSpc];
  odeUserParams.Energy      = new double[nSpc];
  odeUserParams.CvMass      = new double[nSpc];
  odeUserParams.molWt       = new double[nSpc];
  odeUserParams.invMolWt    = new double[nSpc];
  odeUserParams.fwdROP      = new double[nStep];
  odeUserParams.createRate  = new double[nSpc];
  odeUserParams.destroyRate = new double[nSpc];
  odeUserParams.conc        = new double[nSpc];

  // set constant parameters in arrays
  mech->getMolWtSpc(&odeUserParams.molWt[0]);
  for(j=0; j<nSpc; j++) {
    odeUserParams.invMolWt[j]=1.0/odeUserParams.molWt[j];
  }
  // set up the sparse jacobian matrix structure
  odeUserParams.sparseMtx =
    (Jsparse *)alloc_Jsparse(*mech,
                             zerorkCtrl.precThresh,
                             zerorkCtrl.doILU,
                             zerorkCtrl.doFakeUpdate,
                             zerorkCtrl.precThreshType,
                             zerorkCtrl.partialPivotThresh,
                             zerorkCtrl.permutationType);

  // allocate reduced temperature roots
  odeUserParams.nIdtTemp    = cvodeCtrl.nRoots;
  odeUserParams.redTempRoot = new double[cvodeCtrl.nRoots];
  // assign reduced temperature roots
  for(j=0; j<cvodeCtrl.nRoots; j++) {
    odeUserParams.redTempRoot[j] = (parser.idtTemps()[j]+initTemp)/refTemp;
  }
  // sort temperature roots in ascending order
  insertionSort(cvodeCtrl.nRoots,
                &odeUserParams.redTempRoot[0]);

  // set the ropMultiplier pointer
  odeUserParams.ropMultiplierPtr = &ropMultiplier[0];

  // set the pointer to the cvode memory
  odeUserParams.cvodeMemPtr = cvodeCtrl.cvodeMemPtr;

  // allocate and setup the limiter storage

  odeUserParams.step_limiter = new double[nStep];
  for(j=0; j<nStep; ++j) {
    odeUserParams.step_limiter[j] = -1.0; // negative values not used
                                          // for limiter
  }
  odeUserParams.unimolecular_limiter = NULL;
  odeUserParams.bimolecular_limiter  = NULL;

  if(parser.use_unimolecular_limit() ||
     parser.use_bimolecular_limit()) {

    odeUserParams.use_unimolecular_limit = 0;
    odeUserParams.unimolecular_limit = -1.0; // negative values not used
                                             // in limiter function
    if(parser.use_unimolecular_limit()) {
      odeUserParams.use_unimolecular_limit = 1;
      odeUserParams.unimolecular_limit = parser.unimolecular_limit();
    }

    odeUserParams.use_bimolecular_limit = 0;
    odeUserParams.bimolecular_limit = -1.0; // negative values not used
                                             // in limiter function
    if(parser.use_bimolecular_limit()) {
      odeUserParams.use_bimolecular_limit = 1;
      odeUserParams.bimolecular_limit = parser.bimolecular_limit();
    }
    // allocate and initialize storage
    odeUserParams.unimolecular_limiter = new double[nStep];
    odeUserParams.bimolecular_limiter  = new double[nStep];

    for(j=0; j<nStep; ++j) {
      odeUserParams.unimolecular_limiter[j] = -1.0;
      odeUserParams.bimolecular_limiter[j]  = -1.0;
    }

    // store the unimolecular limit, which is fixed
    if(parser.use_unimolecular_limit()) {
      std::vector<int> unimolecular_step_id;
      int num_limit = buildUnimolecularListOfSteps(mech,
                                                   &unimolecular_step_id);

      for(j=0; j<num_limit; ++j) {
        odeUserParams.unimolecular_limiter[unimolecular_step_id[j]] =
          odeUserParams.unimolecular_limit;
      }
    }
  } // end if parser.use_unimolecular_limit() || parser.use_bimolecular_limit()

  // Setup the parameters for the eigenvalue_file
  odeUserParams.compute_eigenvalues = 0;
  if(parser.computeEigenvalues()) {

    odeUserParams.compute_eigenvalues = 1;
    // TODO: fix ugliness
    int name_length = parser.eigenvalueFile().size()+1;
    odeUserParams.eigenvalue_file = new char[name_length];
    strcpy(odeUserParams.eigenvalue_file,parser.eigenvalueFile().c_str());
  }
  else {
    odeUserParams.eigenvalue_file = NULL;
  }

}

void idtControlParams::freeCVodeUserParams(void)
{
  delete [] odeUserParams.redTempRoot;

  free_Jsparse(odeUserParams.sparseMtx);

  delete [] odeUserParams.netProd;
  delete [] odeUserParams.Energy;
  delete [] odeUserParams.CvMass;
  delete [] odeUserParams.molWt;
  delete [] odeUserParams.invMolWt;
  delete [] odeUserParams.fwdROP;
  delete [] odeUserParams.createRate;
  delete [] odeUserParams.destroyRate;
  delete [] odeUserParams.conc;

  if(odeUserParams.unimolecular_limiter != NULL) {
    delete [] odeUserParams.unimolecular_limiter;
  }
  if(odeUserParams.bimolecular_limiter != NULL) {
    delete [] odeUserParams.bimolecular_limiter;
  }
  if(odeUserParams.step_limiter != NULL) {
    delete [] odeUserParams.step_limiter;
  }
  if(odeUserParams.eigenvalue_file != NULL) {
    delete [] odeUserParams.eigenvalue_file;
  }

}


void idtControlParams::setupIdtState(BasicReactorIFP &parser)
{
  int j,spcIdx;
  typedef map< string, double > comp_t;
  comp_t::const_iterator iter;
  double fuelOxygenBal, oxidOxygenBal;
  double *freshMoleFrac, *freshMassFrac;
  double *exhaustMoleFrac, *exhaustMassFrac;

  initTemp = parser.initTemp();
  initPres = parser.initPres();
  initPhi  = parser.initPhi();
  initEgr  = parser.initEgr();
  refTemp  = parser.refTemp();

  // class allocations
  fuelMoleFrac = new double[nSpc];
  oxidMoleFrac = new double[nSpc];
  initMoleFrac = new double[nSpc];
  initMassFrac = new double[nSpc];

  // local allocations
  freshMoleFrac = new double[nSpc];
  freshMassFrac = new double[nSpc];
  exhaustMoleFrac = new double[nSpc];
  exhaustMassFrac = new double[nSpc];

  for(j=0; j<nSpc; j++) {
    fuelMoleFrac[j]=0.0;
    oxidMoleFrac[j]=0.0;
    initMoleFrac[j]=0.0;
    initMassFrac[j]=0.0;
  }

  // set the init fuel composition using the map interator
  for(iter =  parser.fuelComp().begin();
      iter != parser.fuelComp().end();
      ++iter) {

    spcIdx = mech->getIdxFromName(iter->first.c_str());
    if(spcIdx == -1) { 
       printf("ERROR: fuel species %s: not found in mechanism.\n", iter->first.c_str());
       exit(1);
    }
    fuelMoleFrac[spcIdx]=iter->second;
    //printf("fuel species %s (id = %d): %.18g\n",
    //       iter->first.c_str(),spcIdx,fuelMoleFrac[spcIdx]);
  }

  // set the init oxid composition using the map interator
  for(iter =  parser.oxidizerComp().begin();
      iter != parser.oxidizerComp().end();
      ++iter) {

    spcIdx = mech->getIdxFromName(iter->first.c_str());
    if(spcIdx == -1) { 
       printf("ERROR: oxid species %s: not found in mechanism.\n", iter->first.c_str());
       exit(1);
    }
    oxidMoleFrac[spcIdx]=iter->second;
    //printf("oxid species %s (id = %d): %.18g\n",
    //       iter->first.c_str(),spcIdx,oxidMoleFrac[spcIdx]);
  }
  // normalize the mole fraction arrays to one.
  normalizeMoleFrac(&fuelMoleFrac[0]);
  normalizeMoleFrac(&oxidMoleFrac[0]);

  // compute the initial mole fraction of the fresh fuel/oxidizer mixture
  // at the specified equivalence ratio initPhi.
  fuelOxygenBal=mech->getMolarAtomicOxygenRemainder(&fuelMoleFrac[0]);
  oxidOxygenBal=mech->getMolarAtomicOxygenRemainder(&oxidMoleFrac[0]);

  // check for fuel and oxidizer composition that aren't currently handled,
  // and exit the program if found.
  if(fuelOxygenBal > 0.0) {
    printf("ERROR: fuel composition is oxygen rich.\n");
    printf("       moles of atomic-O per mole of fuel is %.18g\n",
           fuelOxygenBal);
    printf("       Contact the Zero-RK developers to add this feature.\n");
    exit(-1);
  }
  if(oxidOxygenBal <= 0.0) {
    printf("ERROR: oxidizer composition contains no oxygen.\n");
    printf("       moles of atomic-O per mole of oxidizer is %.18g\n",
           oxidOxygenBal);
    printf("       Contact the Zero-RK developers to add this feature.\n");
    exit(-1);
  }
  //printf("F/O moles of atomic-O: %.18g/%.18g\n",fuelOxygenBal,oxidOxygenBal);
  for(j=0; j<nSpc; j++) {
    freshMoleFrac[j] = initPhi*fuelMoleFrac[j]
      -fuelOxygenBal/oxidOxygenBal*oxidMoleFrac[j];
  }
  normalizeMoleFrac(&freshMoleFrac[0]);

  // convert the fresh composition to mass fractions
  mech->getYfromX(freshMoleFrac,freshMassFrac);

  // convert the initial composition of the fresh mixture and exhaust
  if(initEgr > 0.0) {
    // get the ideal exhaust molar composition
    mech->getMolarIdealExhaust(freshMoleFrac,exhaustMoleFrac);

    // convert the exhaust mole fractions to mass fractions
    mech->getYfromX(exhaustMoleFrac, exhaustMassFrac);

    // create the initial mass fraction of the blended intake composition
    // here egr represents the fraction by mass the ideal exhaust composition
    // is in the intake composition
    for(j=0; j<nSpc; j++) {
      initMassFrac[j] = (1.0-initEgr)*freshMassFrac[j]
                       +initEgr*exhaustMassFrac[j];
    }
    // convert the initial intake composition to mole fraction
    mech->getXfromY(initMassFrac,initMoleFrac);
  }
  else {
    // zero EGR - copy freshMoleFrac and freshMassFrac
    for(j=0; j<nSpc; j++) {
      initMoleFrac[j]=freshMoleFrac[j];
      initMassFrac[j]=freshMassFrac[j];
    }
  }
  if( parser.pollute_zero_mole_fractions()) {
    for(j=0; j<nSpc; ++j) {
      if(initMoleFrac[j]==0.0) {
         initMoleFrac[j] = parser.absTol()*0.1;
      }
    }
    normalizeMoleFrac(&initMoleFrac[0]);
    mech->getYfromX(initMoleFrac,initMassFrac);
  }

  // set initial density
  initDens=mech->getDensityFromTPY(initTemp, initPres, &initMassFrac[0]);

  delete [] freshMoleFrac;
  delete [] freshMassFrac;
  delete [] exhaustMoleFrac;
  delete [] exhaustMassFrac;
}

void idtControlParams::setupZeroRKParams(BasicReactorIFP &parser)
{
  zerorkCtrl.doFakeUpdate       = parser.doFakeUpdate();
  zerorkCtrl.doILU              = parser.doILU();
  zerorkCtrl.strictSamePattern  = parser.strictSamePattern();

  zerorkCtrl.precThreshType     = parser.precThreshType();
  zerorkCtrl.permutationType    = parser.permutationType();

  zerorkCtrl.precThresh         = parser.precThresh();
  zerorkCtrl.partialPivotThresh = parser.partialPivotThresh();
  zerorkCtrl.permThresh         = parser.permThresh();
}


void idtControlParams::normalizeMoleFrac(double moleFrac[])
{
  int j;
  double sum=0.0;

  for(j=0; j<nSpc; j++) {
    sum+=moleFrac[j];
  }
  if(sum == 0.0) {
    printf("WARNING: In idtControlParams::normalizeMoleFrac(...),\n");
    printf("         can not normalize mole fractions that sum to zero.\n");
  }
  else {
    sum=1.0/sum;
    for(j=0; j<nSpc; j++) {
      moleFrac[j]*=sum;
    }
  }
}
