#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.

#include "utility_funcs.h"

#include "ode_funcs.h"
#include "idtSolvers.h"

using zerork::getHighResolutionTime;

// Computes the ignition delay times for an array of temperature jumps
// for the current ROP multipliers stored in idtCtrl->ropMultiplier[:].
// The temperature jumps are stored in reduced form in the ODE user parameters
// at idtCtrl->odeUserParams.redTempRoot, which is defined by
// odeUserParams.redTempRoot[j] = (parser.idtTemps()[j]+initTemp)/refTemp
// in void idtControlParams::setupCVodeUserParams(AFactorIFP &parser).
//
// The integration ends when either the integration time exceeds 
// parser.maxTime(), or the last root is found with 
// parser.stopAfterLastIdtTemp() set to true. The latter assumes
// the roots are previously sorted in ascending order and that temperature
// increase monotonically.  In situations where temperature is not
// monotonic, multiple root crossings may exists.  However, only the last
// root crossing time is recorded.  More specifically, the last root
// crossing time before the maximum simulation time, CVode failure or the
// last root is found.
//
// results[0] = IDT for first (smallest) temperature jump in parser.idtTemps()
//        [:]
// results[idtCtrl->num_idt_temperatures_]
//            = time for max of first species given in parser.trackSpeciesMax()
// results[idtCtrl->cvodeCtrl->num_roots_]
//            = maximum of first species given in parser.trackSpeciesMax()

int solveIdtSimple(idtControlParams *idtCtrl,
                   double *results,
                   double *solveTime)
{
  //printf("# DEBUG: enter solveIdtSimple(...)\n"); fflush(stdout);
  int j;
  int flag, flagr;
  int errorFlagCtr=0;
  int *rootsFound;
  double startTime=getHighResolutionTime();
  double tcurr = 0.0;
  double tmax  = idtCtrl->cvodeCtrl.maxTime;
  double current_heat_release_rate[2];
  double max_heat_release_rate[2];
  double time_at_max_heat_release_rate[2];
  const int num_track_species_max = idtCtrl->num_track_species_max_;
  const int num_results = idtCtrl->num_results_;
  const int num_roots = idtCtrl->cvodeCtrl.nRoots;
  std::vector<double> species_max;
  std::vector<double> time_at_species_max;
  double *state = NV_DATA_S(idtCtrl->systemState);
  N_Vector chemical_heat_release_rate;
  chemical_heat_release_rate = N_VNew_Serial(2);
  //printf("# DEBUG: solveIdtSimple(...), after variables\n"); fflush(stdout);
  
   
  rootsFound = new int[num_roots];

  max_heat_release_rate[0] = -INFINITY;
  max_heat_release_rate[1] = -INFINITY;
                     
  idtCtrl->resetIdtSim();           // reset CVode

  species_max.assign(num_track_species_max, 0.0);
  time_at_species_max.assign(num_track_species_max, 0.0);
  for(j=0; j<num_track_species_max; ++j) {
    int track_id = idtCtrl->track_species_max_id_[j];
    if(state[track_id] > species_max[j]) {
      species_max[j] = state[track_id];
      time_at_species_max[j] = tcurr;
    }
  }
  // clear the results array
  for(j=0; j<num_results; ++j) {
    results[j] = INFINITY;
  }

  while(tcurr < tmax) {
    flag = CVode(idtCtrl->cvodeCtrl.cvodeMemPtr,
                 tmax,
                 idtCtrl->systemState,
                 &tcurr,
                 CV_ONE_STEP);
    // record heat release ratesand max species if no error (flag < 0) was 
    // found
    if(flag >= 0) {
      int flag_hrr;
      flag_hrr = ChemicalHeatReleaseRate(tcurr,
                                         idtCtrl->systemState,
                                         chemical_heat_release_rate,
                                         &idtCtrl->odeUserParams);

      if(flag_hrr != 0 && idtCtrl->print_level_ > 0) {
        printf(
          "# WARNING: error flag = %d computing ChemicalHeatReleaseRate(...)\n",
          flag_hrr);
        fflush(stdout);
      } else {
        current_heat_release_rate[0] = NV_Ith_S(chemical_heat_release_rate,0);
        current_heat_release_rate[1] = NV_Ith_S(chemical_heat_release_rate,1);

        if(current_heat_release_rate[0] > max_heat_release_rate[0]) {
          max_heat_release_rate[0] = current_heat_release_rate[0];
          time_at_max_heat_release_rate[0] = tcurr;
        }
        if(current_heat_release_rate[1] > max_heat_release_rate[1]) {
          max_heat_release_rate[1] = current_heat_release_rate[1];
          time_at_max_heat_release_rate[1] = tcurr;
        }

      }
      for(j=0; j<num_track_species_max; ++j) {
        int track_id = idtCtrl->track_species_max_id_[j];
        if(state[track_id] > species_max[j]) {
          species_max[j] = state[track_id];
          time_at_species_max[j] = tcurr;
        }
      }
    }

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
          N_VDestroy_Serial(chemical_heat_release_rate);
          return flag;
	}
      }
      else {
        printf("WARNING: CVode() returned error flag=%d\n",flag);
        printf("         at t=%.18g.  Error flag return counter = %d\n",
               tcurr,errorFlagCtr);
        printf("         Stopping solveIdt().\n");
        print_state(idtCtrl->nSpc+1,idtCtrl->systemState);
        N_VDestroy_Serial(chemical_heat_release_rate);
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
        N_VDestroy_Serial(chemical_heat_release_rate);
        return flagr;
      }
      for(j=0; j<num_roots; j++) {
        //if(rootsFound[j] != 0) {
        // ignore species max root for now
        if(rootsFound[j] != 0 && j<idtCtrl->num_idt_temperatures_) {
          results[j] = tcurr;
        }
      }
      if(rootsFound[idtCtrl->num_idt_temperatures_-1] != 0 &&
         idtCtrl->stop_after_last_idt_temp_) {
        // found the last temperature root, so now advance the time past 
        // the stopping time
        tcurr=2.0*tmax;
      }
    }
  }
  // record the IDT time roots from the pointwise maximum
  for(j=0; j<num_track_species_max; ++j) {
    results[j+idtCtrl->num_idt_temperatures_] = time_at_species_max[j];
  }
  // record results values at the last step
  for(j=0; j<num_track_species_max; ++j) {
    results[j+num_roots] = species_max[j];
  }
  results[num_roots+num_track_species_max]   = time_at_max_heat_release_rate[0];
  results[num_roots+num_track_species_max+1] = max_heat_release_rate[0];
  results[num_roots+num_track_species_max+2] = time_at_max_heat_release_rate[1];
  results[num_roots+num_track_species_max+3] = max_heat_release_rate[1];

  // for(j=0; j<num_results; ++j) {
  //   printf("# DEBUG: results[%d] = %22.15e\n", j, results[j]);
  // }

  // if(idtCtrl->print_level_ > 0) {
  //   printf("# Max heat release rate (-E_f(298)*dy/dt): %14.7e at t = %14.7e [s]\n",max_heat_release_rate[0],time_at_max_heat_release_rate[0]);
  //   printf("# Max heat release rate        (-E*dy/dt): %14.7e at t = %14.7e [s]\n",max_heat_release_rate[1],time_at_max_heat_release_rate[1]);

  //   for(j=0; j<num_track_species_max; ++j) {
  //     printf("# Max mass fraction %16s: %14.7e at t = %14.7e [s]\n",
  //            idtCtrl->mech->getSpeciesName(idtCtrl->track_species_max_id_[j]),
  //            species_max[j],
  //            time_at_species_max[j]);
  //   }
    
  // }

  delete [] rootsFound;
  N_VDestroy_Serial(chemical_heat_release_rate);
  *solveTime=getHighResolutionTime()-startTime;
  //printf("# DEBUG: Leaving solveIdtSimple(...)\n");
  return 0;
}

int solveIdtOriginal(idtControlParams *idtCtrl,
                     double *results,
                     double *solveTime)
{
  double startTime = getHighResolutionTime();
  double innerTime;
  int retFlag;
  idtCtrl->clearAllROPMultiplier(); // set all ROP multipliers to one
  retFlag=solveIdtSimple(idtCtrl,
                         results,
                         &innerTime);
  (*solveTime)=getHighResolutionTime()-startTime;
  return retFlag;
}

int solveIdtPerturbRxn(const int rxnId,
                       idtControlParams *idtCtrl,
                       double *results,
                       double *solveTime)
{
  double startTime = getHighResolutionTime();
  double innerTime;
  int retFlag;
  idtCtrl->clearAllROPMultiplier(); // set all ROP multipliers to one

  idtCtrl->setROPMultiplierOfRxn(rxnId,false);
  retFlag=solveIdtSimple(idtCtrl,
                         results,
                         &innerTime);
  idtCtrl->unsetROPMultiplierOfRxn(rxnId);

  if(idtCtrl->doBothDir) {
    // also calculate the idt by dividing the multiplier for the reaction
    idtCtrl->setROPMultiplierOfRxn(rxnId,true);
    retFlag=solveIdtSimple(idtCtrl,
                           &results[idtCtrl->num_results_],
                           &innerTime);
    idtCtrl->unsetROPMultiplierOfRxn(rxnId);
  } 

  (*solveTime)=getHighResolutionTime()-startTime;
  return retFlag;
}
