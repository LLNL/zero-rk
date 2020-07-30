#ifndef IDTCONTROLPARAMS_H
#define IDTCONTROLPARAMS_H

#include <nvector/nvector_serial.h>

#include "zerork/mechanism.h"
#if defined SUNDIALS3
#include <sunlinsol/sunlinsol_dense.h>
#elif defined SUNDIALS4
#include <sunlinsol/sunlinsol_dense.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif


#include "cv_param_sparse.h"
#include "GSA_AFactorIFP.h"

typedef struct {
 
  bool doFakeUpdate; 
  bool doILU;
  bool strictSamePattern;

  int precThreshType;
  int permutationType;

  double precThresh;
  double partialPivotThresh;
  double permThresh;

} ZerorkControlParams;

typedef struct {
  
  int nDims;
  int nRoots;
  int maxSteps;
  int krylovDim;
  int maxNumNonLinIters;

  double relTol;
  double absTol;
  double maxDtInternal;
  double cvEpsLin;
  double cvNlConvCoeff;
  double maxTime;

  void *cvodeMemPtr;
#ifdef SUNDIALS3
  SUNLinearSolver LS;
#elif SUNDIALS4
  SUNLinearSolver LS;
  SUNNonlinearSolver NLS;
#endif

} cvodeControlParams;


class idtControlParams {
 public:
  idtControlParams(char *inpFile, int printLevel);
  ~idtControlParams(void);

  void resetIdtSim(void);
  void clearAllROPMultiplier(void);
  //void setROPMultiplierOfRxn(const int rxnId, const bool doDivision);
  void setROPMultiplierOfRxn(const int rxnId, const double afactor_mult);
  void setROPMultipliers(const double afactor_mult[]);
  void unsetROPMultiplierOfRxn(const int rxnId);
  int nSpc;
  int nRxn;
  int nStep;

  double initTemp;
  double initPres;
  double initPhi;
  double initEgr;

  double initDens;
  double refTemp;

  //bool doBothDir;
  //double AFactorMultiplier;
  double *ropMultiplier;
  N_Vector systemState;

  string mechFile;
  string thermFile;
  string mechLogFile;
  string outFile;
  string checkFile;
  string gsaMatrixFile;

  double *fuelMoleFrac;
  double *oxidMoleFrac;
  double *initMoleFrac;
  double *initMassFrac;

  zerork::mechanism *mech;
  cvodeControlParams cvodeCtrl;
  ZerorkControlParams zerork_ctrl;
  cv_param odeUserParams;


 private:
  void setupZerorkParams(GSA_AFactorIFP &parser);
  void setupIdtState(GSA_AFactorIFP &parser);
  void setupCVode(GSA_AFactorIFP &parser);
  void setupCVodeUserParams(GSA_AFactorIFP &parser);
  void freeCVodeUserParams(void);
  void normalizeMoleFrac(double moleFrac[]);
 

};




#endif
