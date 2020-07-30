#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "utility_funcs.h"

#include "perturbAFactor_common.h"

#include "zerork/constants_api.h"

void checkCommandLine(int inpArgc, char *inpArgv[])
{
  if(inpArgc != 2) {
    printf("ERROR: incorrect command line usage\n");
    printf("       use instead %s <input file>\n",
           inpArgv[0]);
    exit(-1);
  }
  if(access(inpArgv[1],R_OK) != 0) {
    printf("ERROR: can not open input file %s for read access\n",
           inpArgv[1]);
    exit(-2);
  }
}


void getHeaderInfo(int inpArgc,
                   char *inpArgv[],
                   idtControlParams *ctrl,
                   bool doFullReport,
                   string &header)
{
  string str;
  string speciesName;
  header.clear();
  header = "# Command Line:";
  for(int j=0; j<inpArgc; j++) {
    str = inpArgv[j];
    header += " " + str;
  }

  {
    zerork::ZeroRKConstants zrk_constants;
    header += "\n#-------------------------------------------------------------------------------\n";
    header += "# Zero-RK git commit information:\n";
    header += "#    commit id            : " + string(zrk_constants.GetCommitId())        + "\n";
    header += "#    commit time-stamp    : " + string(zrk_constants.GetCommitTimestamp()) + "\n";
    header += "#    branch name          : " + string(zrk_constants.GetBranchName())      + "\n";
  }


  header += "#-------------------------------------------------------------------------------\n";
  header += "# Mechanism Info:\n";
  header += "#    mechanism       filename : " + ctrl->mechFile    + "\n";
  header += "#    thermo          filename : " + ctrl->thermFile   + "\n";
  header += "#    parser log      filename : " + ctrl->mechLogFile + "\n";
  header += "#    A-factor matrix filename : " + ctrl->gsaMatrixFile + "\n";
  header += "#    number of reactions      : " + intToStr(ctrl->nRxn,"%d")  + "\n";
  header += "#    number of steps          : " + intToStr(ctrl->nStep,"%d") + "\n";

  header += "#-------------------------------------------------------------------------------\n";
  header += "# Initial Conditions:\n";
  header += "#    temperature         [K] : "
          + dblToStr(ctrl->initTemp,"%9.4f") + "\n";
  header += "#    pressure           [Pa] : "
    + dblToStr(ctrl->initPres,"%14.7e") + "\n";
  header += "#    equivalence ratio   [-] : "
          + dblToStr(ctrl->initPhi,"%9.7f")  + "\n";
  header += "#    EGR fraction (mass) [-] : "
          + dblToStr(ctrl->initEgr,"%9.7f")  + "\n";
  header += "#-------------------------------------------------------------------------------\n";
  header += "# Initial Composition: (mass fraction, mole fraction)\n";
  for(int j=0; j<ctrl->nSpc; j++) {
    if(ctrl->initMoleFrac[j] > 0.0) {
      speciesName = ctrl->mech->getSpeciesName(j);
      header += "#    species " + intToStr(j+1,"%4d") + " "
	      + speciesName + ": "
	      + dblToStr(ctrl->initMassFrac[j],"%16.7e") + " "
              + dblToStr(ctrl->initMoleFrac[j],"%16.7e") + "\n";
    }
  }
  header += "#-------------------------------------------------------------------------------\n";
  if(doFullReport) {
    header += "# CVode Control Parameters:\n";
    header += "#    relative tolerance         [-] : "
            + dblToStr(ctrl->cvodeCtrl.relTol,"%14.7e") + "\n";
    header += "#    absolute tolerance         [-] : "
            + dblToStr(ctrl->cvodeCtrl.absTol,"%14.7e") + "\n";
    header += "#    max internal dt            [s] : "
            + dblToStr(ctrl->cvodeCtrl.maxDtInternal,"%14.7e") + "\n";
    header += "#    max ode time               [s] : "
            + dblToStr(ctrl->cvodeCtrl.maxTime,"%14.7e") + "\n";
    header += "#    linear solver rel. tol.    [-] : "
            + dblToStr(ctrl->cvodeCtrl.cvEpsLin,"%14.7e") + "\n";
    header += "#    nonlinear solver rel. tol. [-] : "
            + dblToStr(ctrl->cvodeCtrl.cvNlConvCoeff,"%14.7e") + "\n";
    header += "#    max internal steps         [#] : "
            + intToStr(ctrl->cvodeCtrl.maxSteps,"%d") + "\n";
    header += "#    max nonlinear iterations   [#] : "
            + intToStr(ctrl->cvodeCtrl.maxNumNonLinIters,"%d") + "\n";
    header += "#    max Krylov subspace dim.   [#] : "
            + intToStr(ctrl->cvodeCtrl.krylovDim,"%d") + "\n";
    header += "#    number of IDT temperatures [#] : "
            + intToStr(ctrl->cvodeCtrl.nRoots,"%d") + "\n";
    for(int j=0; j<ctrl->cvodeCtrl.nRoots; j++) {
      double temp=ctrl->odeUserParams.redTempRoot[j]*ctrl->refTemp;
      header += "#        T(" + intToStr(j+1,"%2d") + ")   [K] : ";
      header += dblToStr(temp, "%9.4f") + "\n";
    }
    header += "#-------------------------------------------------------------------------------\n";
    header += "# ZeroRK Control Parameters:\n";
    header += "#    preconditioner threshold [-]: "
      + dblToStr(ctrl->zerork_ctrl.precThresh,"%16.7e") + "\n";
    header += "#-------------------------------------------------------------------------------\n";
  } // end doFullReport
}
void getColHeaderInfo(idtControlParams *ctrl,
                      string &header)
{
  header.clear();
  header = "# rxn id    soln time     rem time";
  for(int j=0; j<ctrl->cvodeCtrl.nRoots; j++) {
    header += "   A*mult IDT [s]";
  }
  header += "\n";
  header += "#                 [s]          [s]";
  for(int j=0; j<ctrl->cvodeCtrl.nRoots; j++) {
    double temp=ctrl->odeUserParams.redTempRoot[j]*ctrl->refTemp;
    header += dblToStr(temp, "%13.4f") + " [K]";
  }
  header += "\n";
}

void getColHeaderInfo_IdtTask(idtControlParams *ctrl,
                              string &header)
{
  header.clear();
  header  = "# Ignitional delay times (IDT) based on temperature jump dT\n";
  header += "# Note run id = 0 is the unperturbed IDT calculation\n";
  header += "# -----------------------------------------------------------------------------\n";
  header += "#    run       wall-clock";
  for(int j=0; j<ctrl->cvodeCtrl.nRoots; j++) {
    header += "       IDT [s] at";
  }
  header += "\n";
  header += "#     id         time [s]";
  for(int j=0; j<ctrl->cvodeCtrl.nRoots; j++) {
    double temp=ctrl->odeUserParams.redTempRoot[j]*ctrl->refTemp;
    header += "   dT =" + dblToStr(temp, "%6.1f") + " [K]";
  }
  header += "\n";

}
