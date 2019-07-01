#ifndef CV_PARAM_SPARSE_H
#define CV_PARAM_SPARSE_H

#include "zerork/mechanism.h"

#include "slu_ddefs.h"
#ifdef CVIDT_USE_METIS
#include "metis.h"
#endif




typedef struct
{
  int concIdx;   // or denominator or column index
  int rxnIdx;    // index in the cantera reaction list
  int sparseIdx; // index location in the sparse matrix list
  int sparseIdxCSR; // index location in the sparse matrix list (CSR format)
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
  int nReactors;

  int nSize;
  int nNonZero;
  // compressed column storage
  double *mtxData;  // sparse jacobian data [length nNonZero]
  int *mtxRowIdx;   // row address [length nNonZero]
  int *mtxColSum;   // running tally of elements in each column [length nSize+1]
  int *mtxColIdx;   // CSR format
  int *mtxRowSum;   // CSR format

  // special data access
  int *diagIdx;     // diagonal elements [length nsize]
  int *lastRowIdx;  // last row [length nsize]

  JsparseTermList *termList;

  // sparse solver specific data
  int *isFirstFactor;  // One flag for each reactor

  double offDiagThreshold;

  int permutationType;
  double permThresh;
  double maxGammaChangeOrder;
  bool strictSamePattern;
  int *numPermReUses;
  bool ILU;
  bool fakeUpdate;
  int threshType;

  //char equed[1];
  char* equed;         //nReactors in length
  SuperMatrix* Mslu;
  SuperMatrix* Lslu;
  SuperMatrix* Uslu;
  SuperMatrix* Bslu;
  SuperMatrix* Xslu;
  int *rowPermutation;
  int *colPermutation;
  int *colElimTree; // etree

  int *reduceNNZ;      //nReactors in length
  double *reduceData;
  int *reduceRowIdx;
  int *reduceColSum;
  int *reduceColIdx; //CSR format
  int *reduceRowSum; //CSR format
  int* transform_idx; //CSR format
  int* LUnnz;          //nReactors in length
  double* fillFactor;  //nReactors in length

  double *vecRHS;
  double *vecSoln;
  double *Rvec;
  double *Cvec;
  superlu_options_t* optionSLU;
  SuperLUStat_t* statSLU;
  mem_usage_t*  mem_usage;

  //TODO Testing error metric multipliers for GPU balancing
  double* reactorErrMult;

#ifdef CVIDT_USE_METIS
  //Metis column permutation
  int metis_options[METIS_NOPTIONS];
#endif

#ifdef ZERORK_USE_CUSOLVERRF
  bool cusolverRf;
  int CUSOLVERRF_FACTOR_ALG;
  int CUSOLVERRF_SOLVE_ALG;
  int isFirstFactor_cusolverRf;
  SuperMatrix Mslu_cusolverRf;
  SuperMatrix Bslu_cusolverRf;
  SuperMatrix Xslu_cusolverRf;
  int *rowPermutation_dev;
  int *colPermutation_dev;
  int *reduceRowIdx_dev;
  int *reduceColSum_dev;
  double *reduceData_dev;
  double **reduceDataPtrs;
  double **reduceDataPtrs_dev;
  double ** cusolverRf_rhs_ptrs;
  double ** cusolverRf_rhs_ptrs_dev;
  double * cusolverRf_solve_tmp; //tmp size is bigger for batch routines so need to malloc some space for ourselves
  int *mtxColSum_dev;
  double *mtxData_dev;
  double *diagUnit_dev;
  double *diagUnitCSR_dev;
  int *lastRowIdx_dev;
  // cusolverRf handle data
  cusolverRfHandle_t cusolverRfHandle;
#endif
#ifdef ZERORK_USE_CUSOLVERSP
#ifndef ZERORK_USE_CUSOLVERRF
#error "Must have cusolverRf for cusolverSp"
#endif
  bool cusolverSp;
  int isFirstFactor_cusolverSp;
  double * cusolverSp_solve_tmp;
  int* cusolverSp_perm_map;
  int* csrPermRowSum;
  int* csrPermColIdx;
  // cusolverRf handle data
  cusolverSpHandle_t cusolverSpHandle;
  csrqrInfo_t csrqr_info;
  cusparseMatDescr_t mat_descr;
#endif


  bool abort_matrix_eval;
} Jsparse;

JsparseTermList * alloc_JsparseTermList(const int inpFD, const int inpFC);
void free_JsparseTermList(JsparseTermList *w);

//Jsparse * alloc_Jsparse(const int inpSize, const int inpNNZ, const int inpFD,
//			const int inpFC, const int inpRD, const int inpRC);
Jsparse * alloc_Jsparse(const zerork::mechanism &mechInp, int nReactors);

void free_Jsparse(Jsparse *w);
#ifdef ZERORK_USE_CUSOLVERRF
void Jsparse_setup_cusolverRf(Jsparse *w);
void initialize_cusolverRf(Jsparse *w);
#endif
#ifdef ZERORK_USE_CUSOLVERSP
void Jsparse_setup_cusolverSp(Jsparse *w);
void initialize_cusolverSp(Jsparse *w);
#endif
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
  // reactor type
  bool constPress;

  // df/dt constants
  zerork::mechanism *mech; //static_cast when doing multireactor calls
  int nReactors,currReactor;  // used when mixing single reactor with multi-reactor code
  int nReactorsMax,nReactorsMin;
  int nMatrixReactors;
  int nMatrixThreads; //Ignore OMP_NUM_THREADS for now.
                      // mostly for validation/comparison of GPU vs CPU
  int nSpc;
  double Tref;       // [K] reference temperature
  double Tinit;      // [K] need for root finding on host for single
  double deltaTign;  // [K] delta T used to calculate ignition delay
  double *Press;       // [N/m^2] reactor pressure (used for CP reactors)
  double *invDens;    // [m^3/kg] inverse of the const volume reactor density
  double *molWt;     // [kg/kmol] array of the molecular weights

  // df/dt storage
  double *meanCvMass; // [J/(kg-K)] mixture specific heat at constant volume
  double *dTemp_dt;   // [K/s] temperature time derivative
  double *Energy;    // [J/kg] dimensional internal energy of each species
  double *CvMass;
  double *netProd;
  double *dpdt;

  double *createRate;
  double *destroyRate;
  double *conc;
  double *systemState_host; //Used in to hold data for printing and CPU ops.
  double *systemDeriv_host;
  double *tmp1_host;
  double *tmp2_host;
  double *tmp3_host;
  double *systemDeriv_dev_ptr;
  double *systemState_dev_ptr;
  double *tmp1_dev_ptr;
  double *tmp2_dev_ptr;
  double *tmp3_dev_ptr;
  double *Jac_host;

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

  //TODO:
  int gpu_id;


  // Cvode initialization parameters
  int maxsteps, maxord;
  double rtol, atol, maxdt_internal, nlconvcoef, epslin;
  int abstol_dens; //0 massfrac abstol, 1 molar density abstol

  // use for recreating the divided difference scaling
  void *cvodeMemPtr;
  void *cvodeMemPtr_mr;
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

  //GPU logistic func params
  double logiA, logiK, logiQ, logiM, logiB, logirNu;

  //How much output to include
  int verbosity;

  int rank; //set to 0 if no MPI
  int nranks; //set to 1 if no MPI
#ifdef USE_MPI
  double *rank_weights;
#endif
} cv_param;


#endif
