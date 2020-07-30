#ifndef IDTCONTROLPARAMS_H
#define IDTCONTROLPARAMS_H

#include <vector>

#include <nvector/nvector_serial.h>
#if defined SUNDIALS3
#include <sunlinsol/sunlinsol_dense.h>
#elif defined SUNDIALS4
#include <sunlinsol/sunlinsol_dense.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif

#include "zerork/mechanism.h"

#include "cv_param_sparse.h"
#include "AFactorIFP.h"

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
  // key integers, must be initialized as soon as possible in constructor
  int nSpc;
  int nRxn;
  int nStep;
  int print_level_;
  bool stop_after_last_idt_temp_;
  int num_idt_temperatures_;
  int num_track_species_max_;
  int num_results_;
  std::vector<int> track_species_max_id_;

  double initTemp;
  double initPres;
  double initPhi;
  double initEgr;

  double initDens;
  double refTemp;

  bool doBothDir;
  double AFactorMultiplier;
  double *ropMultiplier;
  N_Vector systemState;

  string mechFile;
  string thermFile;
  string mechLogFile;
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
  void setupZeroRKParams(AFactorIFP &parser);
  void setupIdtState(AFactorIFP &parser);
  void setupCVode(AFactorIFP &parser);
  void setupCVodeUserParams(AFactorIFP &parser);
  void freeCVodeUserParams(void);
  void normalizeMoleFrac(double moleFrac[]);
 

};




#endif
