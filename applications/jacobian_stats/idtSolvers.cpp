#include <math.h>

#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros

#include "ode_funcs.h"
#include "utility_funcs.h"
#include "jacobian_stats.h"

#include "idtSolvers.h"

using zerork::getHighResolutionTime;

// Computes the ignition delay times for an array of temperature jumps
// for the current ROP multipliers stored in idtCtrl->ropMultiplier[:].
// The temperature jumps are stored in reduced form in the ODE user parameters
// at idtCtrl->odeUserParams.redTempRoot, which is defined by
// odeUserParams.redTempRoot[j] = (parser.idtTemps()[j]+initTemp)/refTemp
// in void idtControlParams::setupCVodeUserParams(AFactorIFP &parser).
//
// The integration ends when either the last root is found, this assumes
// the roots are previously sorted in ascending order and that temperature
// increase monotonically.  In situations where temperature is not
// monotonic, multiple root crossings may exists.  However, only the last
// root crossing time is recorded.  More specifically, the last root
// crossing time before the maximum simulation time, CVode failure or the
// last root is found.
int solveIdtSimple(idtControlParams *idtCtrl,
                   double idt[],
                   double *solveTime)
{
  int j;
  int flag, flagr;
  int errorFlagCtr=0;
  int *rootsFound;
  double startTime=getHighResolutionTime();
  double tcurr = 0.0;
  double tmax  = idtCtrl->cvodeCtrl.maxTime;
   
  rootsFound = new int[idtCtrl->cvodeCtrl.nRoots];
  for(j=0; j<idtCtrl->cvodeCtrl.nRoots; j++) {
    idt[j]=INFINITY;
  }

  idtCtrl->resetIdtSim();           // reset CVode

  while(tcurr < tmax) {
    flag = CVode(idtCtrl->cvodeCtrl.cvodeMemPtr,
                 tmax,
                 idtCtrl->systemState,
                 &tcurr,
                 CV_NORMAL);
    // check if error was found
    if(check_flag(&flag, "CVODE ERROR", 1)) {

      errorFlagCtr++;
      if(errorFlagCtr == 1) {
        printf("WARNING: CVode() returned error flag=%d\n",flag);
        printf("         attempting CVodeReInit at t=%.18g [s]\n", tcurr);
        printf("         with preconditioner threshold reset to zero.\n");
	
        change_JsparseThresh(idtCtrl->odeUserParams.sparseMtx,
                             0.0);
        flagr = CVodeReInit(idtCtrl->cvodeCtrl.cvodeMemPtr,
                           tcurr,
                           idtCtrl->systemState);
        if(check_flag(&flagr, "CVodeReInit", 1)) {
          printf("WARNING: CVodeReInit() failed with flag=%d\n",flagr);
          print_state(idtCtrl->nSpc+1,idtCtrl->systemState);
          delete [] rootsFound;
          *solveTime=getHighResolutionTime()-startTime;
          return flag;
	}
      }
      else {
        printf("WARNING: CVode() returned error flag=%d\n",flag);
        printf("         at t=%.18g.  Error flag return counter = %d\n",
               tcurr,errorFlagCtr);
        printf("         Stopping solveIdt().\n");
        print_state(idtCtrl->nSpc+1,idtCtrl->systemState);
        delete [] rootsFound;
        *solveTime=getHighResolutionTime()-startTime;
        return flag;
      }
    }
    else if(flag == CV_ROOT_RETURN) {
      // Scan the roots found and overwrite idt with the current time.
      // This means that if the root is found repeatedly (i.e., the
      // temperature is oscillating) only the time of the last root is
      // recorded. 
      flagr = CVodeGetRootInfo(idtCtrl->cvodeCtrl.cvodeMemPtr,
                               rootsFound);
      if(check_flag(&flagr, "CVodeGetRootInfo", 1)) {
        printf("WARNING: CVodeGetRootInfo returned error flag=%d\n",flagr);
        printf("          at t=%.18g. Stopping solveIdt().\n",tcurr);
        delete [] rootsFound;
        *solveTime=getHighResolutionTime()-startTime;
        return flagr;
      }
      for(j=0; j<idtCtrl->cvodeCtrl.nRoots; j++) {
        if(rootsFound[j] != 0) {
          idt[j] = tcurr;
        }
      }
      if(rootsFound[idtCtrl->cvodeCtrl.nRoots-1] != 0) {
        // found the last root advance the time past the stopping time
        tcurr=2.0*tmax;
      }
    }
  }

  delete [] rootsFound;
  *solveTime=getHighResolutionTime()-startTime;
  return 0;
}

int solveIdtOriginal(idtControlParams *idtCtrl,
                     double idt[],
                     double *solveTime)
{
  double startTime = getHighResolutionTime();
  double innerTime;
  int retFlag;
  idtCtrl->clearAllROPMultiplier(); // set all ROP multipliers to one
  retFlag=solveIdtSimple(idtCtrl,
                         &idt[0],
                         &innerTime);
  (*solveTime)=getHighResolutionTime()-startTime;
  return retFlag;
}

int solveIdtPerturbRxn(const int rxnId,
                       idtControlParams *idtCtrl,
                       double idt[],
                       double *solveTime)
{
  double startTime = getHighResolutionTime();
  double innerTime;
  int retFlag;
  idtCtrl->clearAllROPMultiplier(); // set all ROP multipliers to one

  idtCtrl->setROPMultiplierOfRxn(rxnId,false);
  retFlag=solveIdtSimple(idtCtrl,
                         &idt[0],
                         &innerTime);
  idtCtrl->unsetROPMultiplierOfRxn(rxnId);

  if(idtCtrl->doBothDir) {
    // also calculate the idt by dividing the multiplier for the reaction
    idtCtrl->setROPMultiplierOfRxn(rxnId,true);
    retFlag=solveIdtSimple(idtCtrl,
                           &idt[idtCtrl->cvodeCtrl.nRoots],
                           &innerTime);
    idtCtrl->unsetROPMultiplierOfRxn(rxnId);
  } 

  (*solveTime)=getHighResolutionTime()-startTime;
  return retFlag;
}


int marchIdtStep(idtControlParams *idtCtrl,
                 double t_start,
                 double *t_stop,
                 double *solveTime)
{
  double startTime = getHighResolutionTime();
  double t_current=t_start;
  int flag = CVode(idtCtrl->cvodeCtrl.cvodeMemPtr,
                   (*t_stop),
                   idtCtrl->systemState,
                   &t_current,
                   CV_NORMAL);
  (*t_stop) = t_current;
  (*solveTime) = getHighResolutionTime() - startTime;
  return flag;
}

int marchIdtOneStep(idtControlParams *idtCtrl,
                    double t_start,
                    double *t_stop,
                    double *solveTime)
{
  double startTime = getHighResolutionTime();
  double t_current=t_start;
  int flag = CVode(idtCtrl->cvodeCtrl.cvodeMemPtr,
                   (*t_stop),
                   idtCtrl->systemState,
                   &t_current,
                   CV_ONE_STEP);
  (*t_stop) = t_current;
  (*solveTime) = getHighResolutionTime() - startTime;
  return flag;
}

void initializeIdtStep(idtControlParams *idtCtrl,
                      double idt[])
{
  for(int j=0; j<idtCtrl->cvodeCtrl.nRoots; ++j) {
    idt[j]=INFINITY;
  }

  // Reset the sundials' NVector for systemState, re-initialized cvode
  // to time zero, and set all the solver counters to zero
  idtCtrl->resetIdtSim();
  // set all rate of progress multipliers to 1.0
  idtCtrl->clearAllROPMultiplier();
}


void checkAllFlagsFromIdtStep(idtControlParams *idtCtrl,
                              double t_current,
                              int cvode_flag,
                              double idt[],
                              FILE *check_file)
{
  if(cvode_flag == CV_ROOT_RETURN) {
    // Scan the rootsFound and overwrite idt with the current time.
    // This means that if the root is found repeatedly (i.e., the
    // temperature is oscillating) only the time of the last root is
    // recorded.
    std::vector<int> rootsFound;
    rootsFound.assign(idtCtrl->cvodeCtrl.nRoots,0);
    int flagr = CVodeGetRootInfo(idtCtrl->cvodeCtrl.cvodeMemPtr,
                                 rootsFound.data());
    if(check_flag(&flagr, "CVodeGetRootInfo", 1)) {
       printf("ERROR: CVodeGetRootInfo returned error flag=%d\n",flagr);
       printf("       at t=%.18g. Exiting now!\n",t_current);
       exit(-1);
     }
     for(int j=0; j<idtCtrl->cvodeCtrl.nRoots; ++j) {
       if(rootsFound[j] != 0) {
         idt[j] = t_current;
       }
     }
  }

  if(cvode_flag != CV_SUCCESS && cvode_flag != CV_ROOT_RETURN) {
    // negative values of CVode(...) return flag are errors, while positive
    // values are warnings.
    if(cvode_flag < 0) {
      ++idtCtrl->num_cvode_errors;
    } else {
      ++idtCtrl->num_cvode_warnings;
    }
    writeCVodeWarningsErrors(idtCtrl,
                             t_current,
                             cvode_flag,
                             check_file);  
  }
}

int checkNextIdtStep(idtControlParams *idtCtrl,
                     double t_current,
                     int verbosity)
{
  const int num_species = idtCtrl->mech->getNumSpecies();
  const int num_steps   = idtCtrl->mech->getNumSteps();
  int error_code=0;
  std::vector<double> net_species_rate,step_rate;
  std::string step_description;  
  CVReactorState state;
  getState(idtCtrl,&state);
  getDerivativesFromState(idtCtrl,&state,net_species_rate,step_rate);

  // scan the current state vector for NaNs and inf
  for(int j=0; j<num_species; ++j) {
    if(isfinite(NV_Ith_S(idtCtrl->systemState,j)) == 0) {
      // state vector is either NaN or inf
      printf("Species %d %s isfinite(y) == 0 (false) with y = %lf\n",
             j,
             idtCtrl->mech->getSpeciesName(j),
             NV_Ith_S(idtCtrl->systemState,j));
      fflush(stdout);
      error_code = -1;
    }
  }
  // scan the net reaction rates of the species for NaNs and inf
  for(int j=0; j<num_species; ++j) {
    if(isfinite(net_species_rate[j]) == 0) {
      // state vector is either NaN or inf
      printf("Species %d %s net rate isfinite(ydot) == 0 (false) with ydot = %lf\n",
             j,
             idtCtrl->mech->getSpeciesName(j),
             net_species_rate[j]);
      fflush(stdout);
      error_code += -2;
    }
  }  
  // scan the step rates for NaNs and inf
  for(int j=0; j<num_steps; ++j) {
    if(isfinite(step_rate[j]) == 0) {
      int rxn_id = idtCtrl->mech->getRxnIdxOfStep(j);
      idtCtrl->mech->getReactionNameDirOfStep(j,&step_description);
      printf("Step %d rate isfinite(rop) == 0 (false) with rop = %lf\n",
             j,
             step_rate[j]);
      printf("  from reaction %d: %s\n", rxn_id, step_description.c_str());
      fflush(stdout);
    }
  }

  return error_code;
}



void getState(idtControlParams *idtCtrl,
              CVReactorState *state)
{
  state->mole_fraction.assign(idtCtrl->mech->getNumSpecies(),0.0);
  state->mass_fraction.assign(idtCtrl->mech->getNumSpecies(),0.0);
  state->concentration.assign(idtCtrl->mech->getNumSpecies(),0.0);

  // assign the mass fractions
  for(int j=0; j<idtCtrl->mech->getNumSpecies(); ++j) {
    state->mass_fraction[j] = NV_Ith_S(idtCtrl->systemState,j);
  }

  // convert the mass fractions to mole fractions
  idtCtrl->mech->getXfromY(state->mass_fraction.data(),
                           state->mole_fraction.data());
  // convert the mass fractions to concentration
  idtCtrl->mech->getCfromVY(1.0/idtCtrl->initDens,
                            state->mass_fraction.data(),
                            state->concentration.data());

  state->density = idtCtrl->initDens; // constant volume reactor
  state->temperature = idtCtrl->refTemp*
                       NV_Ith_S(idtCtrl->systemState,
                       idtCtrl->mech->getNumSpecies());
  state->pressure = idtCtrl->mech->getPressureFromTVY(state->temperature,
                                                  1.0/idtCtrl->initDens,
						  state->mass_fraction.data());
  state->molecular_weight = 
    idtCtrl->mech->getMolWtMixFromY(state->mass_fraction.data());
}

// Prerequisite - getState called to initialize the CVReactorState
void getDerivativesFromState(idtControlParams *idtCtrl,
                             CVReactorState *state,
                             std::vector<double> &net_species_rate,
                             std::vector<double> &step_rate)
{
  const int num_species = idtCtrl->mech->getNumSpecies();
  const int num_steps   = idtCtrl->mech->getNumSteps();
  std::vector<double> creation_rate, destruction_rate;

  net_species_rate.clear();
  step_rate.clear();
  
  net_species_rate.assign(num_species,0.0);
  step_rate.assign(num_steps,0.0);
  creation_rate.assign(num_species,0.0);
  destruction_rate.assign(num_species,0.0);

  idtCtrl->mech->getReactionRatesLimiter(state->temperature,
                                  &state->concentration[0],
				  idtCtrl->odeUserParams.step_limiter,
                                  &net_species_rate[0],
                                  &creation_rate[0],
                                  &destruction_rate[0],
                                  &step_rate[0]);
}

void writeCVodeWarningsErrors(idtControlParams *idtCtrl,
                              double t_current,
                              int cvode_flag,
                              FILE *check_file)
{
  writeCVodeReport(idtCtrl,
                   cvode_flag,
                   check_file);
  writeThermoStateReport(idtCtrl,
                         t_current,
                         check_file);

  writeOdeStateReport(idtCtrl,
                      t_current,
                      check_file);
  fprintfDashedLine(check_file);
  writeROPDerivativeReport(idtCtrl,
                           t_current,
                           -1.0, // minimum value to report
                           check_file);
  fprintfDashedLine(check_file);
  writeROPTempDerivativeReport(idtCtrl,
                               t_current,
                               -1.0,
                               check_file);
   // state vector, rhs, normalized rhs?
  // error sorting
  // jacobian element sorting
  // rhs discontinuity (fwd step, rev time step)?
  fprintf(check_file,"\n\n");
}
void writeCVodeReport(idtControlParams *idtCtrl,
                    int cvode_flag,
                    FILE *fptr)
{
  int method_order;
  double time_step;
  char *flag_name;

  fprintfDashedLine(fptr);
  if(cvode_flag < 0) {
    fprintf(fptr,
            "# CVode Error %d\n",
            idtCtrl->num_cvode_errors);
  } else {
    fprintf(fptr,
            "# CVode Warning %d\n",
            idtCtrl->num_cvode_warnings);
  }
  // CVodeGetReturnFlagName(...) is allocating a character array with malloc
  // for its return value, which means the user needs to free() it to avoid
  // a memory leak and to run valgrind clean.
  flag_name = CVodeGetReturnFlagName(cvode_flag);
  fprintf(fptr,
          "# CVode(...) return flag = %d (%s)\n",
          cvode_flag,
          flag_name);
  free(flag_name); // free character array from CVodeGetReturnFlagName(...)
  
  // last order, last step refers to the last internal timestep taken
  CVodeGetLastOrder(idtCtrl->cvodeCtrl.cvodeMemPtr,&method_order);
  CVodeGetLastStep(idtCtrl->cvodeCtrl.cvodeMemPtr,&time_step);
  fprintf(fptr,
          "#   Method order of the last internal timestep         : %d\n",
          method_order);
  fprintf(fptr,
          "#   Step size of the last internal timestep        [s] : %.18g \n",
          time_step);
  // current order, current step refers to the next internal timestep that
  // is going to be attempted
  CVodeGetCurrentOrder(idtCtrl->cvodeCtrl.cvodeMemPtr,&method_order);
  CVodeGetCurrentStep(idtCtrl->cvodeCtrl.cvodeMemPtr,&time_step);
  fprintf(fptr,
          "#   Method order to try on the next internal timestep  : %d\n",
          method_order);
  fprintf(fptr,
          "#   Step size to try on the next internal timestep [s] : %.18g \n",
          time_step);
}
void writeThermoStateReport(idtControlParams *idtCtrl,
                            double t_current,
                            FILE *fptr)
{
  CVReactorState current_state;
  double total_conc;

  getState(idtCtrl,&current_state);

  total_conc = current_state.pressure/
    (idtCtrl->mech->getGasConstant()*current_state.temperature);

  fprintfDashedLine(fptr);
  fprintf(fptr,
          "# Thermodynamic State at Time %.18g [s]\n",
          t_current);
  fprintf(fptr,
          "#   temperature                [K] : %.18g\n",
          current_state.temperature);
  fprintf(fptr,
          "#   pressure                  [Pa] : %.18g\n",
          current_state.pressure);
  fprintf(fptr,
          "#   density               [kg/m^3] : %.18g\n",
          current_state.density);
  fprintf(fptr,
          "#   total concentration [kmol/m^3] : %.18g\n",
          total_conc);
  fprintf(fptr,
          "#   molecular_weight     [kg/kmol] : %.18g\n",
          current_state.molecular_weight);
  fprintf(fptr,
          "#                                                               [kmol/m^3]\n");
  fprintf(fptr,
          "#   index,        species,   mole fraction,   mass fraction,   concentration\n");
  
  for(int j=0; j<idtCtrl->mech->getNumSpecies(); ++j) {
    fprintf(fptr,
            "%8d %16s %16.7e %16.7e %16.7e\n",
            j,
            idtCtrl->mech->getSpeciesName(j),
            current_state.mole_fraction[j],
            current_state.mass_fraction[j],
            current_state.concentration[j]);
  }           
}

void writeOdeStateReport(idtControlParams *idtCtrl,
                         double t_current,
                         FILE *fptr)
{
  int temperature_id = idtCtrl->mech->getNumSpecies();
  N_Vector state_derivative, error_weight, error_estimate;
  state_derivative = N_VNew_Serial(idtCtrl->cvodeCtrl.nDims);
  error_weight     = N_VNew_Serial(idtCtrl->cvodeCtrl.nDims);
  error_estimate   = N_VNew_Serial(idtCtrl->cvodeCtrl.nDims);

  CVodeGetErrWeights(idtCtrl->cvodeCtrl.cvodeMemPtr,
                     error_weight);
  CVodeGetEstLocalErrors(idtCtrl->cvodeCtrl.cvodeMemPtr,
                         error_estimate);
  const_vol_wsr_limiter(t_current,
                        idtCtrl->systemState,
                        state_derivative,
                        &(idtCtrl->odeUserParams));
                        
  fprintfDashedLine(fptr);
  fprintf(fptr,
          "# CVode ODE State at Time %.18g [s]\n",
          t_current);
  fprintf(fptr,
          "# Reference temperature T_ref = %.18g [K]\n",
          idtCtrl->refTemp);
  fprintf(fptr,
          "# The failure score is the absolute value of the product of the local error \n# estimate and the error weights.  The largest value of this product may      \n# indicate the component that is causing recent error test failures according \n# to the CVode User's Guide.\n#\n");
  fprintf(fptr,
          "#                                                                    error      failure\n");
  fprintf(fptr,
          "#   index               description   state vector     derivative   estimate     score\n");

  for(int j=0; j<idtCtrl->mech->getNumSpecies(); ++j) {
    fprintf(fptr,
            "%8d %16s-mass frac %14.7e %14.7e %10.3e %10.3e\n",
            j,
            idtCtrl->mech->getSpeciesName(j),
            NV_Ith_S(idtCtrl->systemState,j),
            NV_Ith_S(state_derivative,j),
            NV_Ith_S(error_estimate,j),
            fabs(NV_Ith_S(error_estimate,j)*NV_Ith_S(error_weight,j)));
  }  
  // print out the temperature information
  fprintf(fptr,
          "%8d                    T/T_ref %14.7e %14.7e %10.3e %10.3e\n",
          temperature_id,
          NV_Ith_S(idtCtrl->systemState,temperature_id),
          NV_Ith_S(state_derivative,temperature_id),
          NV_Ith_S(error_estimate,temperature_id),
          fabs(NV_Ith_S(error_estimate,temperature_id)*
               NV_Ith_S(error_weight,temperature_id)));
         

  N_VDestroy_Serial(state_derivative);
  N_VDestroy_Serial(error_weight);
  N_VDestroy_Serial(error_estimate);
  
}

void printfDashedLine()
{
  printf("# ----------------------------------------------------------------------------\n");
} 
void fprintfDashedLine(FILE *fptr)
{
  fprintf(fptr,"# ----------------------------------------------------------------------------\n");
} 

double GetMinMoleFraction(const CVReactorState &state,
                          int *min_species_id)
{
  double min_value = 1.0e+300;
  *min_species_id  = -1;

  for(size_t j=0; j<state.mole_fraction.size(); ++j) {
    if(state.mole_fraction[j] < min_value) {
      *min_species_id = (int)j;
      min_value = state.mole_fraction[j];
    }
  }
  return min_value;
}

double GetNegativeMoleFractionSum(const CVReactorState &state)
{
  double sum = 0.0;
  for(size_t j=0; j<state.mole_fraction.size(); ++j) {
    if(state.mole_fraction[j] < 0.0) {
      sum += state.mole_fraction[j];
    }
  }
  return sum;
}
