#ifndef IDTCONTROLPARAMS_H
#define IDTCONTROLPARAMS_H

#include <stdio.h>
#include <nvector/nvector_serial.h>
#if defined SUNDIALS3
#include <sunlinsol/sunlinsol_dense.h>
#elif defined SUNDIALS4
#include <sunlinsol/sunlinsol_dense.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif


#include "zerork/mechanism.h"

#include "cv_param_sparse.h"
#include "BasicReactorIFP.h"

typedef struct {

  bool doFakeUpdate;
  bool doILU;
  bool strictSamePattern;

  int precThreshType;
  int permutationType;

  double precThresh;
  double partialPivotThresh;
  double permThresh;

} zerorkControlParams;

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
  void setROPMultiplierOfRxn(const int rxnId, const bool doDivision);
  void unsetROPMultiplierOfRxn(const int rxnId);
  int nSpc;
  int nRxn;
  int nStep;

  int num_cvode_warnings;
  int num_cvode_errors;

  double initTemp;
  double initPres;
  double initPhi;
  double initEgr;

  double initDens;
  double refTemp;

  double printDeltaTime;

  int num_scans;
  double scan_temperature_ref;
  double min_scan_temperature;
  double max_scan_temperature;

  bool doBothDir;
  double AFactorMultiplier;
  double *ropMultiplier;
  N_Vector systemState;

  string mechFile;
  string thermFile;
  string mechLogFile;
  string mechStatFile;
  string jacobianStatFile;
  string jacobianRawFile;
  string integratorProbFile;
  string outFile;

  double *fuelMoleFrac;
  double *oxidMoleFrac;
  double *initMoleFrac;
  double *initMassFrac;

  zerork::mechanism *mech;
  cvodeControlParams cvodeCtrl;
  zerorkControlParams zerorkCtrl;
  cv_param odeUserParams;


 private:
  void setupZeroRKParams(BasicReactorIFP &parser);
  void setupIdtState(BasicReactorIFP &parser);
  void setupCVode(BasicReactorIFP &parser);
  void setupCVodeUserParams(BasicReactorIFP &parser);
  void freeCVodeUserParams(void);
  void normalizeMoleFrac(double moleFrac[]);


};




#endif
