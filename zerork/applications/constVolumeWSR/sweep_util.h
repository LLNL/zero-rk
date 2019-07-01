#ifndef SWEEP_UTIL_H
#define SWEEP_UTIL_H

#include "zerork/mechanism.h"

const int MAX_SPEC_NAME_LEN = 32;
const int MAX_SPEC_MOLE_LEN = 32; 
const int MAX_FILE_NAME_LEN = 1024;

double readFileLine_first_double(FILE *fptr);
int readFileLine_first_int(FILE *fptr);
char * readFileLine_first_string(FILE *fptr);
int readFileLine_first_comp(FILE *fptr, zerork::mechanism &mechInp,
			    double moleFrac[]);

const char permutationTypeName[3][27]=
  {"MMD_AT_PLUS_A",
   "Metis_NodeND",
   "External (from perm_c.txt)"};
const char precondThresholdName[5][24]=
  {"absolute",
   "column-norm",
   "row-norm",
   "min(row,col) norm",
   "fixed fraction of terms"};


class idt_sweep_params
{
 public:
  idt_sweep_params(char *inputFileName);
  void setInitialComp(const double phi, const double egr);
  void printInitialMoleComp() const;
  void incrementRunId();
  
  double getInitTemp()  const {return initTemp[tempId];}
  double getInitPres()  const {return initPres[presId];}
  double getPhi()       const {return phi[phiId];}
  double getEgr()       const {return egr[egrId];}
  double getThresh()    const {return thresh[threshId];}
  int    getKrylovDim() const {return krylovDim[krylovId];}
  int    getThreshId()  const {return threshId;}
  double getThreshFromId(const int id)  const {return thresh[id];}
  int    getNumThreshRuns() const {return nThreshRuns;}
  int    getRunId()     const {return runId;}
  int    getRunTotal()  const {return runTotal;}
  double getDensity()   const;
  void   getInitMoleFrac(double x[]) const;
  void   getInitMassFrac(double y[]) const;
  int    getNumTrackSpecies() const {return nTrackSpc;}
  int    getTrackSpeciesId(const int id) const {return trackSpcIdx[id];}

  int    getNumFuelSpecies() const {return nFuelSpc;}
  int    getNumOxidSpecies() const {return nOxidSpc;}
  int    getFuelSpeciesId(const int id) const {return fuelSpcIdx[id];}
  int    getOxidSpeciesId(const int id) const {return oxidSpcIdx[id];}
  double getFuelMoleFrac(const int id) const {return fuelMoleFrac[id];}
  double getOxidMoleFrac(const int id) const {return oxidMoleFrac[id];}

  int    getNumSpecies() const {return nSpc;}
  int    getNumSteps()   const {return gasMech->getNumSteps();}
  double getRefTemp()    const {return Tref;}
  double getDeltaTign()  const {return deltaTign;}
  double getRTol()       const {return rtol;}
  double getATol()       const {return atol;}
  int    getMaxInternalSteps() const {return maxInternalSteps;}
  double getMaxInternalDt()    const {return maxInternalDt;}
  double getStopTime()         const {return stopTime;}
  double getPrintTime()        const {return printTime;}
  char * getIdtFileName()      const {return idtFileName;}
  char * getTHistFileName()    const {return thistFileName;}
  char * getMechFileName()     const {return mechFileName;}
  char * getThermFileName()    const {return thermFileName;}
  char * getLogFileName()      const {return logFileName;}

  bool getILU()         const {return doILU;}
  int getPermutationType()         const {return permutationType;}
  bool getUpdate()        const {return fakeUpdate;}
  int getThreshType()        const {return threshType;}
  double getPartialPivotThresh()        const {return partialPivotThresh;}
  double getEpsLin()        const {return epsLin;}
  double getNlConvCoeff()        const {return nlConvCoeff;}
  bool oneStep()         const {return oneStepMode;}
  //bool printAll()        const {return printAllSteps;}
  int continueAfterIDT()  const {return isContinueAfterIDT;}

  zerork::mechanism * getMechPtr() const {return gasMech;}

  int getMaxPrimaryCvodeFails() const {return maxCvodeFails1;}
  int getMaxSecondaryCvodeFails() const {return maxCvodeFails2;}
  double getSafetyThreshold() const {return safetyThresh;}

  double getFuelOxyBalance() const {return fuelOxygenBal;}
  double getOxidOxyBalance() const {return oxidOxygenBal;}
  double getMoleOxidStoicPerFuel() const {return moleOxidStoicPerFuel;}

  ~idt_sweep_params();

 private:
  // data directly read from input file
  int nFuelSpc,nOxidSpc;
  int isConstPres;
  int isContinueAfterIDT;
  double stopTime,printTime,maxInternalDt;
  int maxInternalSteps;
  double rtol,atol,deltaTign,Tref;
  bool fakeUpdate,doILU,oneStepMode;
  //bool printAllSteps;
  int permutationType;
  int threshType;
  double partialPivotThresh,epsLin,nlConvCoeff;

  int nTempRuns,nPresRuns,nPhiRuns,nEgrRuns,nThreshRuns,nKrylovRuns;
  // ----------------------------------
  int nSpc, nTrackSpc;
  double fuelOxygenBal,oxidOxygenBal,moleOxidStoicPerFuel;

  char *mechFileName,*thermFileName,*logFileName,*idtFileName,*thistFileName;

  zerork::mechanism *gasMech;
  
  double *fuelMoleFrac;
  double *oxidMoleFrac;
  double *freshMoleFrac,   *freshMassFrac;
  double *exhaustMoleFrac, *exhaustMassFrac;
  double *initMoleFrac,    *initMassFrac;

  int *trackSpcIdx,*fuelSpcIdx,*oxidSpcIdx;

  int tempId,presId,phiId,egrId,threshId,krylovId,runId;
  int runTotal;
  double *initTemp,*initPres, *phi,*egr,*thresh;
  int *krylovDim;

  int maxCvodeFails1,maxCvodeFails2;
  double safetyThresh;

  double normalize(const int N, double v[]);
};

#endif
