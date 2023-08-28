#include <stdlib.h>
#include <stdio.h>
#ifdef _WIN32
#include <io.h>
#define F_OK 00
#define R_OK 04
#define access _access
#else
#include "unistd.h"
#endif

#include <vector>

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
  header += "#    mechanism  file name : " + ctrl->mechFile    + "\n";
  header += "#    thermo     file name : " + ctrl->thermFile   + "\n";
  header += "#    parser log file name : " + ctrl->mechLogFile + "\n";
  header += "#    number of species    : " + intToStr(ctrl->nSpc,"%d")  + "\n";
  header += "#    number of reactions  : " + intToStr(ctrl->nRxn,"%d")  + "\n";
  header += "#    number of steps      : " + intToStr(ctrl->nStep,"%d") + "\n";
  header += "#    A-factor multiplier  : "
          + dblToStr(ctrl->AFactorMultiplier,"%6.4f") + "\n";
  if(ctrl->doBothDir) {
    header += "#    perturb both dirs    : yes (mult & div)\n";
  } else {
    header += "#    perturb both dirs    : no  (mult only)\n";
  }

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
            + intToStr(ctrl->num_idt_temperatures_,"%d") + "\n";
    for(int j=0; j<ctrl->num_idt_temperatures_; j++) {
      double temp=ctrl->odeUserParams.redTempRoot[j]*ctrl->refTemp;
      header += "#        T(" + intToStr(j+1,"%2d") + ")   [K] : ";
      header += dblToStr(temp, "%9.4f") + "\n";
    }
    header += "#    number of species to track max [#] : "
            + intToStr(ctrl->num_track_species_max_,"%d") + " (species index starts at zero)\n";
    for(int j=0; j<ctrl->num_track_species_max_; j++) {
      int track_id = ctrl->track_species_max_id_[j];
      header += "#        y(" + intToStr(track_id,"%d") + ") : ";
      header += std::string(ctrl->mech->getSpeciesName(track_id)) + "\n";
    }
    header += "#-------------------------------------------------------------------------------\n";
    header += "# Zero-RK Control Parameters:\n";
    header += "#    preconditioner threshold [-]: "
      + dblToStr(ctrl->zerorkCtrl.precThresh,"%16.7e") + "\n"; 
    header += "#-------------------------------------------------------------------------------\n";
  } // end doFullReport
}
void getColHeaderInfo(idtControlParams *ctrl,
                      string &header)
{
  const int COLUMN_PAD = 21;
  int offset;
  //const int num_results = ctrl->num_results_;
  std::vector<std::string> max_species_names;

  for(int j=0; j<ctrl->num_track_species_max_; ++j) {
    int track_id = ctrl->track_species_max_id_[j];
    max_species_names.push_back(ctrl->mech->getSpeciesName(track_id));
  }
  // Create enumerated header list to facilitate plotting
  header.clear();
  header  = "# Column  1: reaction index (0 is the unperturbed solution)\n";
  header += "# Column  2: [s] CPU run time\n";
  header += 
    "# Column  3: [s] estimated remaining run time (interactive cases)\n";

  //--------------------------------------------------------------------------
  // IDT for temperature jump
  offset = 4;
  for(int j=0; j<ctrl->num_idt_temperatures_; j++) {
    header += "# Column " + intToStr(offset,"%2d") + 
      ": [s] A*mult IDT at T = " + 
      dblToStr(ctrl->odeUserParams.redTempRoot[j]*ctrl->refTemp,"%13.4f") + 
      " [K]\n";
    ++offset;
  }
  // IDT for max species
  for(int j=0; j<ctrl->num_track_species_max_; ++j) {
    header += "# Column " + intToStr(offset,"%2d") +
      ": [s] A*mult time at max mass fraction of " + max_species_names[j]+"\n";
    ++offset;
  }
  // max species value
  for(int j=0; j<ctrl->num_track_species_max_; ++j) {
    header += "# Column " + intToStr(offset,"%2d") +
      ": [-] A*mult max mass fraction of " + max_species_names[j]+"\n";
    ++offset;
  }
  // max heat release rate
  header += "# Column " + intToStr(offset,"%2d") +
    ": [s] A*mult time at max heat release rate type 1 (HRR1 uses energy of formation at T = 298.15)\n";
  ++offset;
  header += "# Column " + intToStr(offset,"%2d") +
    ": [J/m^3/s] A*mult max heat release rate type 1 (HRR1 uses energy of formation at T = 298.15)\n";
  ++offset;
  header += "# Column " + intToStr(offset,"%2d") +
    ": [s] A*mult time at max heat release rate type 2 (HRR2 uses internal energy at the current temperature)\n";
  ++offset;
  header += "# Column " + intToStr(offset,"%2d") +
    ": [J/m^3/s] A*mult max heat release rate type 2 (HRR2 uses internal energy at the current temperature)\n";
  ++offset;

  if(ctrl->doBothDir) {
    // IDT for temperature jump
    for(int j=0; j<ctrl->num_idt_temperatures_; j++) {
      header += "# Column " + intToStr(offset,"%2d") + 
        ": [s] A/mult IDT at T = " + 
        dblToStr(ctrl->odeUserParams.redTempRoot[j]*ctrl->refTemp,"%13.4f") + 
        " [K]\n";
      ++offset;
    }
    // IDT for max species
    for(int j=0; j<ctrl->num_track_species_max_; ++j) {
      header += "# Column " + intToStr(offset,"%2d") +
        ": [s] A/mult time at max mass fraction of " + max_species_names[j]+"\n";
      ++offset;
    }
    // max species value
    for(int j=0; j<ctrl->num_track_species_max_; ++j) {
      header += "# Column " + intToStr(offset,"%2d") +
        ": [-] A/mult max mass fraction of " + max_species_names[j]+"\n";
      ++offset;
    }
    // max heat release rate
    header += "# Column " + intToStr(offset,"%2d") +
      ": [s] A/mult time at max heat release rate type 1 (HRR1 uses energy of formation at T = 298.15)\n";
    ++offset;
    header += "# Column " + intToStr(offset,"%2d") +
      ": [J/m^3/s] A/mult max heat release rate type 1 (HRR1 uses energy of formation at T = 298.15)\n";
    ++offset;
    header += "# Column " + intToStr(offset,"%2d") +
      ": [s] A/mult time at max heat release rate type 2 (HRR2 uses internal energy at the current temperature)\n";
    ++offset;
    header += "# Column " + intToStr(offset,"%2d") +
      ": [J/m^3/s] A/mult max heat release rate type 2 (HRR2 uses internal energy at the current temperature)\n";
    ++offset;
  }
  if(ctrl->doBothDir) {
    header += "# Column " + intToStr(offset,"%2d") +
      ": [J/m^3/s] Delta HRR1 (A*mult - A/mult)\n";
    ++offset;
    header += "# Column " + intToStr(offset,"%2d") +
      ": [J/m^3/s] Delta HRR2 (A*mult - A/mult)\n";
    ++offset;
  } else {
    header += "# Column " + intToStr(offset,"%2d") +
      ": [J/m^3/s] Delta HRR1 (A*mult - A(orig))\n";
    ++offset;
    header += "# Column " + intToStr(offset,"%2d") +
      ": [J/m^3/s] Delta HRR2 (A*mult - A(orig))\n";
    ++offset;
  }



  //--------------------------------------------------------------------------
  // Create column by column listing
  // pad the species strings by the column pad
  for(int j=0; j<ctrl->num_track_species_max_; ++j) {
    int num_pad = COLUMN_PAD - max_species_names[j].size(); 
    if(num_pad > 0) {
      max_species_names[j].insert(0,num_pad, ' ');
    }
  }

  header += "# rxn id    soln time     rem time";
  for(int j=0; j<ctrl->num_idt_temperatures_; j++) {
    header += "       A*mult IDT [s]";
  }
  for(int j=0; j<ctrl->num_track_species_max_; j++) {
    header += "     A*mult t@max [s]";
  }
  for(int j=0; j<ctrl->num_track_species_max_; j++) {
    //         123456789012345678901
    header += "    max mass frac [-]";
  }
  //         123456789012345678901
  header += "     A*mult t@max [s]";
  header += "     A*mult [J/m^3/s]";
  header += "     A*mult t@max [s]";
  header += "     A*mult [J/m^3/s]";

  if(ctrl->doBothDir) {
    for(int j=0; j<ctrl->num_idt_temperatures_; j++) {
      header += "       A/mult IDT [s]";
    }
    for(int j=0; j<ctrl->num_track_species_max_; j++) {
      header += "     A/mult t@max [s]";
    }
    for(int j=0; j<ctrl->num_track_species_max_; j++) {
      //         123456789012345678901
      header += "    max mass frac [-]";
    }
    //         123456789012345678901
    header += "     A/mult t@max [s]";
    header += "     A/mult [J/m^3/s]";
    header += "     A/mult t@max [s]";
    header += "     A/mult [J/m^3/s]";
  }
  if(ctrl->doBothDir) {
    //         123456789012345678901
    header += "      A*mult - A/mult";
    header += "      A*mult - A/mult";
  } else {
    //         123456789012345678901
    header += "     A*mult - A(orig)";
    header += "     A*mult - A(orig)";
  }

  header += "\n";
  header += "#                 [s]          [s]";
  for(int j=0; j<ctrl->num_idt_temperatures_; j++) {
    double temp=ctrl->odeUserParams.redTempRoot[j]*ctrl->refTemp;
    header += dblToStr(temp, "    %13.4f") + " [K]";
  }
  for(int j=0; j<ctrl->num_track_species_max_; j++) {
    header += max_species_names[j];
  }
  for(int j=0; j<ctrl->num_track_species_max_; j++) {
    header += max_species_names[j];
  }
  //         123456789012345678901
  header += "             max HRR1";
  header += "             max HRR1";
  header += "             max HRR2";
  header += "             max HRR2";

  if(ctrl->doBothDir) {
    for(int j=0; j<ctrl->num_idt_temperatures_; j++) {
      double temp=ctrl->odeUserParams.redTempRoot[j]*ctrl->refTemp;
      header += dblToStr(temp, "    %13.4f") + " [K]";
    }
    for(int j=0; j<ctrl->num_track_species_max_; j++) {
      header += max_species_names[j];
    }
    for(int j=0; j<ctrl->num_track_species_max_; j++) {
      header += max_species_names[j];
    }
    //         123456789012345678901
    header += "             max HRR1";
    header += "             max HRR1";
    header += "             max HRR2";
    header += "             max HRR2";

  }    
  //         123456789012345678901
  header += " delta HRR1 [J/m^3/s]";
  header += " delta HRR2 [J/m^3/s]";
  header += "\n";
}

void getColHeaderInfo_sensitivity(idtControlParams *ctrl,
                                  string &header)
{
  header.clear();
  header = "# Relative sensitivity: Srel = d(ln(IDT)) / d(ln(k)) = (k/IDT)*d(IDT)/d(k)\n";
  header += "#      approximated by: Srel = ln(IDT_2/IDT_1)/ln(k_2/k_1)\n";
  header += "#\n";
  header += "#      IDT_2 is the ignition delay after multiplying the rate of progress\n";
  header += "#      by the A-factor multiplier.  IDT_1 is the ignition delay of the original\n";
  header += "#      mechanism (unperturbed), or the IDT after dividing the rate of progress\n";
  header += "#      by the A-factor multiplier, depending on if doBothDir in the input file\n";
  header += "#      is set to 'n' or 'y' respectively.\n";
  header += "#------------------------------------------------------------------------------\n";
  header += "# rxn id";
  for(int j=0; j<ctrl->num_idt_temperatures_; j++) {
    header += "      Srel [-] at";
  }
  header += "\n";
  header += "#       ";
  for(int j=0; j<ctrl->num_idt_temperatures_; j++) {
    double temp=ctrl->odeUserParams.redTempRoot[j]*ctrl->refTemp;
    header += dblToStr(temp, "%13.4f") + " [K]";
  }
  header += "\n";   
}

// idtOrig[] length ctrl->cvodeCtrl.nRoots*nRxn
// idtPerturb[] length nRxn*ctrl->cvodeCtrl.nRoots if ctrl->doBothDir == false
//                   2*nRxn*ctrl->cvodeCtrl.nRoots if ctrl->doBothDir == true
// rxnSens[] length nRxn*ctrl->num_idt_temperatures_
void calcRxnSensitivity(const idtControlParams *ctrl,
                        const double idtOrig[],
                        const double idtPerturb[],
                        double rxnSens[],
                        int sortedRxnIdx[])
{
  int nRxn = ctrl->mech->getNumReactions();
  int nIdtTemp = ctrl->num_idt_temperatures_;
  int num_solutions = ctrl->num_results_;
  double currMax,currSens;
  maxRxnSens_t *maxRxnSensList;

  maxRxnSensList = new maxRxnSens_t[nRxn]; 

  for(int j=0; j<nRxn; j++) {
    currMax = 0.0;
    for(int k=0; k<nIdtTemp; k++) {
      if(ctrl->doBothDir) {
        currSens = log(idtPerturb[j*2*num_solutions+k]/
                       idtPerturb[j*2*num_solutions+num_solutions+k])/
	           log(2.0*ctrl->AFactorMultiplier);
      }
      else {
        currSens = log(idtPerturb[j*num_solutions+k]/idtOrig[k])/
	           log(ctrl->AFactorMultiplier);
      }
      if(fabs(currSens) > currMax) {
        currMax = fabs(currSens);
      }
      rxnSens[j*nIdtTemp+k] = currSens;
    }
    maxRxnSensList[j].rxnId = j;
    maxRxnSensList[j].maxRelSens = currMax;
  }

  // sort by the largest relative sensitivity magnitude
  qsort((void *)&maxRxnSensList[0],
        nRxn,
        sizeof(maxRxnSens_t),
        compare_maxRxnSens_t);

  for(int j=0; j<nRxn; j++) {
    sortedRxnIdx[j] = maxRxnSensList[j].rxnId;
  }
  delete [] maxRxnSensList;
}

int compare_maxRxnSens_t(const void *A, const void *B)
{
  maxRxnSens_t *Aptr =(maxRxnSens_t *)A;
  maxRxnSens_t *Bptr =(maxRxnSens_t *)B;
  if(Aptr->maxRelSens < Bptr->maxRelSens) {
    return 1;
  }
  else if (Aptr->maxRelSens > Bptr->maxRelSens) {
    return -1;
  }
  return 0;
}
