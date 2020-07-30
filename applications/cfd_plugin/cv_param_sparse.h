#ifndef CV_PARAM_SPARSE_H
#define CV_PARAM_SPARSE_H

#include "zerork/mechanism.h"

#include "slu_ddefs.h"

#if defined SUNDIALS3
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#elif defined SUNDIALS4
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif

typedef struct
{
  int concIdx;   // or denominator or column index
  int rxnIdx;    // index in the cantera reaction list
  int sparseIdx; // index location in the sparse matrix list
} JsparseIdx;

typedef struct
{
  int nFwdDestroy;  // FD suffix
  int nFwdCreate;   // FC suffix

  JsparseIdx *fwdDestroy;
  JsparseIdx *fwdCreate;

} JsparseTermList;

typedef struct
{
  int nSize;
  int nNonZero;
  // compressed column storage
  double *mtxData;  // sparse jacobian data [length nNonZero]
  int *mtxRowIdx;   // row address [length nNonZero]
  int *mtxColSum;   // running tally of elements in each column [length nSize+1]

  // special data access
  int *diagIdx;     // diagonal elements [length nsize]
  int *lastRowIdx;  // last row [length nsize]

  JsparseTermList *termList;

  int isFirstFactor;

  double offDiagThreshold;

  double permThresh;
  double maxGammaChangeOrder;
  bool strictSamePattern;
  int numPermReUses;
  bool fakeUpdate;

  char equed[1];
  SuperMatrix Mslu;
  SuperMatrix Lslu;
  SuperMatrix Uslu;
  SuperMatrix Bslu;
  SuperMatrix Xslu;
  int *rowPermutation;
  int *colPermutation;
  int *colElimTree; // etree

  int reduceNNZ;
  double *reduceData;
  int *reduceRowIdx;
  int *reduceColSum;
  int LUnnz;
  double fillFactor;

  double *vecRHS;
  double *vecSoln;
  double *Rvec;
  double *Cvec;
  superlu_options_t optionSLU;
  SuperLUStat_t statSLU;
  mem_usage_t mem_usage;
#if SUPERLU_MAJOR_VERSION > 4
  GlobalLU_t Glu;
#endif

  bool abort_matrix_eval;
} Jsparse;

JsparseTermList * alloc_JsparseTermList(const int inpFD, const int inpFC);
void free_JsparseTermList(JsparseTermList *w);

Jsparse * alloc_Jsparse(const zerork::mechanism &mechInp);
void free_Jsparse(Jsparse *w);
void print_Jsparse(Jsparse *w);
double getElement_Jsparse(Jsparse *w,const int rowIdx, const int colIdx);
void calcReaction_Jsparse(Jsparse *w, const double invPosConc[],
			  const double fwdROP[]);
void change_JsparseThresh(Jsparse *w, double newThresh);


void countJacobianTerms(const zerork::mechanism &mechInp, int *nFwdDestroy, int *nFwdCreate);

typedef struct
{
  int nFwdDestroy;  // FD suffix
  int nFwdCreate;   // FC suffix

  int *rxnFD, *rxnFC;
  int *specNumFD, *specNumFC;
  int *specDenFD, *specDenFC;

} Jterm_param;

Jterm_param * alloc_Jterm_param(zerork::mechanism &mechInp);
void free_Jterm_param(Jterm_param *w);
void print_Jterm_param(Jterm_param *w);

typedef double jacreal; // variable type for jacobian processing

typedef struct
{
  // df/dt constants
  zerork::mechanism *mech; //static_cast when doing multireactor calls
  int nSpc;
  double Tref;       // [K] reference temperature
  double Tinit;      // [K] need for root finding on host for single
  double deltaTign;  // [K] delta T used to calculate ignition delay
  double Press;       // [N/m^2] reactor pressure (used for CP reactors)
  double invDens;    // [m^3/kg] inverse of the const volume reactor density
  double *molWt;     // [kg/kmol] array of the molecular weights

  // df/dt storage
  double meanCvMass; // [J/(kg-K)] mixture specific heat at constant volume
  double dTemp_dt;   // [K/s] temperature time derivative
  double *Energy;    // [J/kg] dimensional internal energy of each species
  double *CvMass;
  double *netProd;

  double *createRate;
  double *destroyRate;
  double *conc;

  // Jacobian constants
  double minMassFrac;
  double sqrtUnitRnd;
  double *invMolWt;

  // Jacobian storage
  bool doingJacSetup;
  double *fwdROP;
  //jacreal *Mmtx,*Bvec,*Xvec; // ,*Jsave;
  //int *ipvt;

  // sparse Jacobian storage
  Jsparse *sparseMtx;

  // Cvode initialization parameters
  int maxsteps, maxord;
  double rtol, atol, maxdt_internal, nlconvcoef, epslin;
  int abstol_dens; //0 massfrac abstol, 1 molar density abstol

  void *cvodeMemPtr;
#if defined SUNDIALS3
  SUNMatrix A;
  SUNLinearSolver LS;
#elif defined SUNDIALS4
  SUNMatrix A;
  SUNLinearSolver LS;
  SUNNonlinearSolver NLS;
#endif

  // cvode stats
  long int nsteps;
  long int nfevals;
  long int nlinsetups;
  long int netfails;
  long int nniters;
  long int nncfails;
  long int ngevals;
  int nFunc,nJacSetup,nJacRescale,nJacFactor,nBackSolve,nColPerm;
  double colPermTime,jacSetupTime,precSetupTime,jacFactorTime,
         backsolveTime,funcTime,simTime;
  int lsolver; //which solver to use (saved for GPU initialization)

  long int prevNumErrTestFails;

  //How much output to include
  int verbosity;

  int rank; //set to 0 if no MPI
  int nranks; //set to 1 if no MPI
#ifdef USE_MPI
  double *rank_weights;
#endif
} cv_param;


#endif
