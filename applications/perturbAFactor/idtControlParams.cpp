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


idtControlParams::idtControlParams(char *inpFile, int printLevel)
{
  int j;
  AFactorIFP parser(inpFile);
  const double standard_temperature = 298.15;
  double *enthalpy_of_formation;
  double *standard_entropy;
  double *molecular_mass;

  // -----------------------------------------------------------------------
  // zerork mechanism setup:
  mechFile    = parser.mechFile();
  thermFile   = parser.thermFile();
  mechLogFile = parser.mechLogFile();
  outFile     = parser.outFile();

  print_level_ = printLevel;
  stop_after_last_idt_temp_ = parser.stopAfterLastIdtTemp();
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
  num_idt_temperatures_ = (int)parser.idtTemps().size();

  // set up the track species max array
  track_species_max_id_.clear();
  num_track_species_max_ = (int)parser.trackSpeciesMax().size();
  // must be adjusted at the end to ignore missing species

  for(j=0; j<num_track_species_max_; ++j) {
    int species_id = mech->getIdxFromName(parser.trackSpeciesMax()[j].c_str());
    if(0 <= species_id && species_id < nSpc) {
      track_species_max_id_.push_back(species_id);
    } else if(print_level_ > 0) {
      // species not found
      printf("# WARNING: skipping input from trackSpeciesMax[%d] = %s.\n"
             "#          Not found in the mechanism.\n",
             j,
             parser.trackSpeciesMax()[j].c_str());
      fflush(stdout);
    }
  }
  // update with the actual number of valid species
  num_track_species_max_ = (int)track_species_max_id_.size();

  // the results array will carry the IDT for the user specified temperature
  // jumps the maximum value and its corresponding time found for any species
  // mass fractions set in the input file vector trackSpeciesMax, and
  // the maximum heat release rate and time of the maximum heat release rate
  // computed via two methods
  num_results_ = num_idt_temperatures_ + 2*num_track_species_max_ + 4;

  // get the AFactor perturbation settings.
  doBothDir = parser.doBothDir();
  AFactorMultiplier = parser.AFactorMultiplier();
  if(AFactorMultiplier == 1.0) {
    printf("ERROR: AFactorMultiper can not equal one.\n");
    exit(-1);
  }

  ropMultiplier = new double[nStep];
  clearAllROPMultiplier();  // set ropMultiplier[:] = 1.0

  const double gas_constant = mech->getGasConstant();
  odeUserParams.energy_of_formation_ = new double[nSpc];
  enthalpy_of_formation = new double[nSpc];
  standard_entropy = new double[nSpc];

  mech->getIntEnergy_RT(standard_temperature, odeUserParams.energy_of_formation_);
  mech->getEnthalpy_RT(standard_temperature, enthalpy_of_formation);
  //  G/RT = H/RT - S/R  =>  S/R = H/RT - G/RT
  mech->getNonDimGibbsFromT(standard_temperature, standard_entropy);

  for(j=0; j<nSpc; ++j) {
    // convert to S/R
    standard_entropy[j] = enthalpy_of_formation[j]-standard_entropy[j];
    // convert to J/kmol and J/kmol/K
    odeUserParams.energy_of_formation_[j]  *= gas_constant*standard_temperature;
    enthalpy_of_formation[j] *= gas_constant*standard_temperature;
    standard_entropy[j]      *= gas_constant;
  }
  if(print_level_ > 0) {
    FILE *log_file_ptr = fopen(mechLogFile.c_str(),"a");
    if(log_file_ptr == NULL) {
      printf("# Thermodynamic reference properties at T = %8.3f [K]\n"
             "#                  Energy of   Enthalpy of   Standard\n"
             "#                  Formation    Formation    Entropy\n"
             "#   Species Name   [kJ/mol]     [kJ/mol]     [J/mol/K]\n",
             standard_temperature);

      for(j=0; j<nSpc; ++j) {
        printf("%16s  %9.3f    %9.3f    %9.3f\n",
               mech->getSpeciesName(j),
               odeUserParams.energy_of_formation_[j]*1.0e-6,
               enthalpy_of_formation[j]*1.0e-6,
               standard_entropy[j]*1.0e-3);
      }
      fflush(stdout);
    } else {
      fprintf(log_file_ptr,
              "# Thermodynamic reference properties at T = %8.3f [K]\n"
              "#                  Energy of   Enthalpy of   Standard\n"
              "#                  Formation    Formation    Entropy\n"
              "#   Species Name   [kJ/mol]     [kJ/mol]     [J/mol/K]\n",
              standard_temperature);

      for(j=0; j<nSpc; ++j) {
        fprintf(log_file_ptr,
                "%16s  %9.3f    %9.3f    %9.3f\n",
                mech->getSpeciesName(j),
                odeUserParams.energy_of_formation_[j]*1.0e-6,
                enthalpy_of_formation[j]*1.0e-6,
                standard_entropy[j]*1.0e-3);
      }
      fflush(log_file_ptr);
      fclose(log_file_ptr);
    }
  }
  // convert energy of formation from molar to mass basis
  molecular_mass = standard_entropy; // use allocated space
  mech->getMolWtSpc(molecular_mass);
  for(int j=0; j<nSpc; ++j) {
    odeUserParams.energy_of_formation_[j] *= 1.0/molecular_mass[j]; // [J/kmol/K] to [J/kg/K]
  }


  setupZeroRKParams(parser);
  setupCVode(parser);       // includes allocating N_Vector systemState
  setupIdtState(parser);
  setupCVodeUserParams(parser);
  resetIdtSim(); // sets up the system state

  // clean up temporary arrays
  delete [] standard_entropy;
  delete [] enthalpy_of_formation;
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

void idtControlParams::setupCVode(AFactorIFP &parser)
{
  int flag;

  // set cvode control parameters
  cvodeCtrl.nDims             = nSpc+1;
  cvodeCtrl.nRoots            = num_idt_temperatures_+num_track_species_max_;
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
                 const_vol_wsr_perturb,
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
  // tempRootFunc[0] = (T_idt[0] - T)/T_ref,  want zero-crossing from + to -
  //             [:]
  // tempRootFunc[num_idt_temperatures_] = ydot[track_species_max_id_[0]
  //                                       want zero-crossing from + to -
  //                                       for maxima
  std::vector<int> root_directions;
  for(int j=0; j<cvodeCtrl.nRoots; ++j) {
    root_directions.push_back(-1);
  }

  flag = CVodeRootInit(cvodeCtrl.cvodeMemPtr,
                       cvodeCtrl.nRoots,
                       tempRootFunc);
  if (check_flag(&flag, "CVodeRootInit", 1)) {
    exit(-1);
  }
  /* Set the root directions to only catch + to - zero-crossings */
  flag = CVodeSetRootDirection(cvodeCtrl.cvodeMemPtr,
                               &root_directions[0]);
  if (check_flag(&flag, "CVodeSetRootDirection", 1)) {
    exit(-1);
  }

 /* Disable the rootfinding warnings when problem initially seems to be
    at a root for multiple steps */
  flag = CVodeSetNoInactiveRootWarn(cvodeCtrl.cvodeMemPtr);
  if (check_flag(&flag, "CVodeSetNoInactiveRootWarn", 1)) {
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

void idtControlParams::setupCVodeUserParams(AFactorIFP &parser)
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

  odeUserParams.state_workspace_ = N_VNew_Serial(nSpc+1);

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
  odeUserParams.nIdtTemp    = num_idt_temperatures_;
  odeUserParams.redTempRoot = new double[num_idt_temperatures_];
  // allocate track species max array
  odeUserParams.num_track_species_max_ = num_track_species_max_;
  odeUserParams.track_species_max_id_ = NULL;

  if(num_track_species_max_ > 0) {
    odeUserParams.track_species_max_id_ = new int[num_track_species_max_];
  }

  for(j=0; j<num_track_species_max_; ++j) {
    odeUserParams.track_species_max_id_[j] = track_species_max_id_[j];
  }

  // assign reduced temperature roots
  for(j=0; j<num_idt_temperatures_; j++) {
    odeUserParams.redTempRoot[j] = (parser.idtTemps()[j]+initTemp)/refTemp;
  }
  // sort temperature roots in ascending order
  insertionSort(num_idt_temperatures_,
                &odeUserParams.redTempRoot[0]);

  // set the ropMultiplier pointer
  odeUserParams.ropMultiplierPtr = &ropMultiplier[0];

  // set the pointer to the cvode memory
  odeUserParams.cvodeMemPtr = cvodeCtrl.cvodeMemPtr;
}

void idtControlParams::freeCVodeUserParams(void)
{
  delete [] odeUserParams.redTempRoot;

  free_Jsparse(odeUserParams.sparseMtx);

  delete [] odeUserParams.energy_of_formation_;
  delete [] odeUserParams.netProd;
  delete [] odeUserParams.Energy;
  delete [] odeUserParams.CvMass;
  delete [] odeUserParams.molWt;
  delete [] odeUserParams.invMolWt;
  delete [] odeUserParams.fwdROP;
  delete [] odeUserParams.createRate;
  delete [] odeUserParams.destroyRate;
  delete [] odeUserParams.conc;


  N_VDestroy_Serial(odeUserParams.state_workspace_);

  if(odeUserParams.track_species_max_id_ != NULL) {
    delete [] odeUserParams.track_species_max_id_;
  }

}


void idtControlParams::setupIdtState(AFactorIFP &parser)
{
  int j,spcIdx;
  typedef map< string, double > comp_t;
  comp_t::const_iterator iter;
  double fuelOxygenBal, oxidOxygenBal;
  double *freshMoleFrac, *freshMassFrac;
  double *exhaustMoleFrac, *exhaustMassFrac;

  // TODO: change to vector input for T, P, phi
  if(parser.initTemp().size() > 1 ||
     parser.initPres().size() > 1 ||
     parser.initPhi().size() > 1) {

    if(print_level_ > 0) {
      printf("# WARNING: Support for vector valued initial conditions not implemented.\n"
             "#          only the first values are used for {T_0, p_0, and phi_0}:\n"
             "#              T_0   = %24.17e [K]\n"
             "#              p_0   = %24.17e [Pa]\n"
             "#              phi_0 = %24.17e [-]\n",
             parser.initTemp()[0],
             parser.initPres()[0],
             parser.initPhi()[0]);
      fflush(stdout);
    }
  }

  initTemp = parser.initTemp()[0];
  initPres = parser.initPres()[0];
  initPhi  = parser.initPhi()[0];
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
  if(fuelOxygenBal > 0.0 && initPhi > 0.0) {
    printf("ERROR: fuel composition is oxygen rich.\n");
    printf("       moles of atomic-O per mole of fuel is %.18g\n",
           fuelOxygenBal);
    printf("       Contact the Zero-RK developers to add this feature.\n");
    exit(-1);
  }
  if(oxidOxygenBal <= 0.0 && initPhi > 0.0) {
    printf("ERROR: oxidizer composition contains no oxygen.\n");
    printf("       moles of atomic-O per mole of oxidizer is %.18g\n",
           oxidOxygenBal);
    printf("       Contact the Zero-RK developers to add this feature.\n");
    exit(-1);
  }
  if(initPhi <= 0.0 && initEgr > 0.0) {
    printf("ERROR: can not create an initial composition with\n");
    printf("       equivalence ratio (F/A)/(F/A)_stoic = %.18g\n",
           initPhi);
    printf("       egr mass fraction = %.18g\n",
           initEgr);
    printf("       Contact the Zero-RK developers to add this capability.\n");
    exit(-1);
  }

  if(initPhi > 0.0) {
    for(j=0; j<nSpc; j++) {
      freshMoleFrac[j] = initPhi*fuelMoleFrac[j]
        -fuelOxygenBal/oxidOxygenBal*oxidMoleFrac[j];
    }
    normalizeMoleFrac(&freshMoleFrac[0]);
  }
  else { // equivalence ratio is zero and the composition is only the
         // oxidizer
    printf("WARNING: equivalence ratio (F/A)/(F/A)_stoic  is zero\n");
    printf("         setting the initial composition to the oxidizer\n");
    for(j=0; j<nSpc; j++) {
      freshMoleFrac[j] = oxidMoleFrac[j];
    }
  }

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

  // set initial density
  initDens=mech->getDensityFromTPY(initTemp, initPres, &initMassFrac[0]);

  delete [] freshMoleFrac;
  delete [] freshMassFrac;
  delete [] exhaustMoleFrac;
  delete [] exhaustMassFrac;
}

void idtControlParams::setupZeroRKParams(AFactorIFP &parser)
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
