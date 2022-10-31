#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.

#include "utilities.h"

#include "utility_funcs.h"

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

int solveIdtGsaPerturbRxn(const double afactor_mult[],
                          idtControlParams *idtCtrl,
                          double idt[],
                          double *solveTime)
{
  double startTime = getHighResolutionTime();
  double innerTime;
  int retFlag;

  idtCtrl->setROPMultipliers(&afactor_mult[0]);
  retFlag=solveIdtSimple(idtCtrl,
                         &idt[0],
                         &innerTime);
  (*solveTime)=getHighResolutionTime()-startTime;
  return retFlag;
}

//int solveIdtPerturbRxn(const int rxnId,
//                       idtControlParams *idtCtrl,
//                       double idt[],
//                       double *solveTime)
//{
//  double startTime = getHighResolutionTime();
//  double innerTime;
//  int retFlag;
//  idtCtrl->clearAllROPMultiplier(); // set all ROP multipliers to one
//
//  idtCtrl->setROPMultiplierOfRxn(rxnId,false);
//  retFlag=solveIdtSimple(idtCtrl,
//                         &idt[0],
//                         &innerTime);
//  idtCtrl->unsetROPMultiplierOfRxn(rxnId);
//
//  if(idtCtrl->doBothDir) {
//    // also calculate the idt by dividing the multiplier for the reaction
//    idtCtrl->setROPMultiplierOfRxn(rxnId,true);
//    retFlag=solveIdtSimple(idtCtrl,
//                           &idt[idtCtrl->cvodeCtrl.nRoots],
//                           &innerTime);
//    idtCtrl->unsetROPMultiplierOfRxn(rxnId);
//  } 
//
//  (*solveTime)=getHighResolutionTime()-startTime;
//  return retFlag;
//}
