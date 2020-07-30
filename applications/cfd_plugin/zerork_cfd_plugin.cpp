
//System libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <time.h>
#include <sys/time.h>
#include <algorithm>
#include <unistd.h> //for geteuid
#include <cfloat> //for DBL_MAX
#include <numeric> //for accumulate

#include "zerork_cfd_plugin.h"
#include "cv_param_sparse.h"
#include "ode_funcs.h"
#include "matrix_funcs.h"
#include "utility_funcs.h"

//Sundials
#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.


#if defined SUNDIALS2
#include <cvode/cvode_dense.h>      // prototypes & constants for CVDense
#include <cvode/cvode_spgmr.h>      // prototypes & constants for CVSPGMR
#elif defined SUNDIALS3
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <cvode/cvode_spils.h>
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spgmr.h>
#elif defined SUNDIALS4
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif

#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros


#ifdef USE_MPI
#include "mpi.h"
#endif

#ifdef ZERORK_TRAP_FE
#include <fenv.h> //for fpe trapping
#endif

#define ZERORK_WORK_SHARING

//GLOBAL VAR
namespace {
cv_param systemParam;
FILE* reactor_log_file;
} // anonymous namespace


static FILE* open_reactor_log_file(const char* filename);
static void zerork_update_cvode_stats(cv_param* cvp, void* cvode_mem_ptr, bool multireac);

//Used by LLNL converge udf
extern "C"
void zerork_cfd_plugin_setup_full
(
    int verbosity,
    int maxsteps,
    int maxord,
    double rtol,
    double atol,
    int abstol_dens,
    double maxdt_internal,
    double sparsethresh,
    double nlconvcoef,
    double epslin,
    int lsolver,
    const char* mech_filename,
    const char* therm_filename,
    const char* ck_log_filename,
    const char* reactor_log_filename,
    int* multireac,
    int* gpu_id
)
{
  int j;
  int nSpc,nStep,nState;
  // cvode variables
  int flag;
  // timing data
  double startTime,stopTime;
  double setupCVode,setupCantera,setupJterm;
  zerork::mechanism *gasMech;

#ifdef ZERORK_TRAP_FE
  //set fpe on
  //feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW );
#endif

#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&(systemParam.rank));
  MPI_Comm_size(MPI_COMM_WORLD,&(systemParam.nranks));
#else
  systemParam.rank = 0;
  systemParam.nranks = 1;
#endif

  // open ignition delay file for write and print header
  reactor_log_file = open_reactor_log_file(reactor_log_filename);
  if(systemParam.rank == 0)
  {
      fprintf(reactor_log_file,"#%12s","solve_number");
      fprintf(reactor_log_file," %16s %16s %16s %16s",
                               "reactors_solved","n_cpu","n_gpu",
                               "n_gpu_groups");
      fprintf(reactor_log_file," %16s %16s %16s %16s",
                               "n_steps_avg", "n_steps_avg_cpu",
                               "n_steps_avg_gpu","max_time_cpu");
      fprintf(reactor_log_file," %16s %16s %16s %16s\n",
                              "max_time_gpu","step_time_cpu",
                              "step_time_gpu", "max_time_total");
      fflush(reactor_log_file);
  }

  *gpu_id = -1;
  *multireac = 0;

  int zerork_verbosity = 0;
  if(systemParam.rank==0) zerork_verbosity = 1;
  // Build zerork::mechanism
  gasMech = new zerork::mechanism(mech_filename,therm_filename,
                                  ck_log_filename,zerork_verbosity);

  nSpc=gasMech->getNumSpecies();
  nStep=gasMech->getNumSteps();
  nState=nSpc+1;

  // set up the system parameters
  systemParam.maxsteps=maxsteps;
  systemParam.maxord=maxord;
  systemParam.rtol=rtol;
  systemParam.atol=atol;
  systemParam.abstol_dens = abstol_dens;
  systemParam.maxdt_internal=maxdt_internal;
  systemParam.nlconvcoef=nlconvcoef;
  systemParam.epslin=epslin;
  systemParam.verbosity = verbosity;
  systemParam.lsolver = lsolver;

  systemParam.doingJacSetup = false;
  systemParam.nSpc=nSpc;
  systemParam.minMassFrac=1.0e-30;
  systemParam.sqrtUnitRnd=sqrt(UNIT_ROUNDOFF);
  systemParam.Tref=1.0;
  systemParam.mech=gasMech;

  systemParam.netProd    =(double *)malloc(sizeof(double)*nSpc);
  systemParam.Energy     =(double *)malloc(sizeof(double)*nSpc);
  systemParam.CvMass     =(double *)malloc(sizeof(double)*nSpc);
  systemParam.molWt      =(double *)malloc(sizeof(double)*nSpc);
  systemParam.invMolWt   =(double *)malloc(sizeof(double)*nSpc);
  systemParam.fwdROP     =(double *)malloc(sizeof(double)*nStep);
  systemParam.createRate =(double *)malloc(sizeof(double)*nSpc);
  systemParam.destroyRate=(double *)malloc(sizeof(double)*nSpc);
  systemParam.conc       =(double *)malloc(sizeof(double)*nSpc);

  zerork_reset_stats();

  // set constant parameters
  gasMech->getMolWtSpc(systemParam.molWt);
  for(j=0; j<nSpc; j++)
    {systemParam.invMolWt[j]=1.0/systemParam.molWt[j];}

  systemParam.sparseMtx= (Jsparse *)alloc_Jsparse(*(systemParam.mech));

  systemParam.sparseMtx->permThresh = 0.3;
  systemParam.sparseMtx->maxGammaChangeOrder = 3.0;
  systemParam.sparseMtx->strictSamePattern = 0;

  systemParam.sparseMtx->offDiagThreshold = sparsethresh;
  systemParam.sparseMtx->fakeUpdate = false;
  systemParam.sparseMtx->optionSLU.DiagPivotThresh = 0.0;

  systemParam.sparseMtx->reduceNNZ=0; // won't be set until the first J
  systemParam.sparseMtx->LUnnz = 0; // won't be set until the first J
  systemParam.sparseMtx->fillFactor= 100.; // won't be set until the first J
  systemParam.sparseMtx->isFirstFactor = 1;
  systemParam.sparseMtx->numPermReUses = 0;
  systemParam.prevNumErrTestFails = 0;

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration */
#ifdef SUNDIALS4
  systemParam.cvodeMemPtr = CVodeCreate(CV_BDF);
#else
  systemParam.cvodeMemPtr = CVodeCreate(CV_BDF, CV_NEWTON);
#endif
  if (check_cvode_flag(systemParam.cvodeMemPtr, "CVodeCreate", 0)) exit(-1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  N_Vector y0 = N_VNew_Serial(nState);
  flag=CVodeInit(systemParam.cvodeMemPtr, const_vol_wsr, 0.0, y0);
  if (check_cvode_flag(&flag, "CVodeInit", 1)) exit(-1);

  //N.B.:  We re-set these tolerances at every solve init.
  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  flag = CVodeSStolerances(systemParam.cvodeMemPtr, rtol, atol);
  if (check_cvode_flag(&flag, "CVodeSStolerances", 1)) exit(-1);

  /* Call CVodeRootInit to specify the root function with 1 component */
//  flag = CVodeRootInit(systemParam.cvodeMemPtr, 1, tempRootFunc);
//  //flag = CVodeRootInit(systemParam.cvodeMemPtr, 2, tempRoot2);
//  if (check_cvode_flag(&flag, "CVodeRootInit", 1)) exit(-1);

  /* Set the pointer to user-defined data */
  flag = CVodeSetUserData(systemParam.cvodeMemPtr, &systemParam);
  if(check_cvode_flag(&flag, "CVodeSetUserData", 1)) exit(-1);

  if(lsolver == 0 || lsolver == 1) {
#if defined SUNDIALS2
    flag = CVDense(systemParam.cvodeMemPtr, nState);
    if(check_cvode_flag(&flag, "CVDense", 1)) exit(-1);

    if(lsolver == 0) {
      flag = CVDlsSetDenseJacFn(systemParam.cvodeMemPtr, jac_full_dense);
      if(check_cvode_flag(&flag, "CVDlsDenseJacFn", 1)) exit(-1);
    } else if (lsolver == 1) {
      flag = CVDlsSetDenseJacFn(systemParam.cvodeMemPtr, NULL);
      if(check_cvode_flag(&flag, "CVDlsDenseJacFn", 1)) exit(-1);
    }
#elif defined SUNDIALS3
    /* Create dense SUNMatrix for use in linear solves */
    systemParam.A = SUNDenseMatrix(nState, nState);
    if(check_cvode_flag((void *)systemParam.A, "SUNDenseMatrix", 0)) exit(-1);

    /* Create dense SUNLinearSolver object for use by CVode */
    systemParam.LS = SUNDenseLinearSolver(y0, systemParam.A);
    if(check_cvode_flag((void *)systemParam.LS, "SUNDenseLinearSolver", 0)) exit(-1);

    /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
    flag = CVDlsSetLinearSolver(systemParam.cvodeMemPtr, systemParam.LS, systemParam.A);
    if(check_cvode_flag(&flag, "CVDlsSetLinearSolver", 1)) exit(-1);

    if(lsolver == 0) {
      /* Set the user-supplied Jacobian routine Jac */
      flag = CVDlsSetJacFn(systemParam.cvodeMemPtr, jac_full_dense);
      if(check_cvode_flag(&flag, "CVDlsSetJacFn", 1)) exit(-1);
    } else if (lsolver == 1) {
      flag = CVDlsSetJacFn(systemParam.cvodeMemPtr, NULL);
      if(check_cvode_flag(&flag, "CVDlsSetJacFn", 1)) exit(-1);
    }
#elif defined SUNDIALS4
    /* Create dense SUNMatrix for use in linear solves */
    systemParam.A = SUNDenseMatrix(nState, nState);
    if(check_cvode_flag((void *)systemParam.A, "SUNDenseMatrix", 0)) exit(-1);

    /* Create dense SUNLinearSolver object for use by CVode */
    systemParam.LS = SUNDenseLinearSolver(y0, systemParam.A);
    if(check_cvode_flag((void *)systemParam.LS, "SUNDenseLinearSolver", 0)) exit(-1);

    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
    flag = CVodeSetLinearSolver(systemParam.cvodeMemPtr, systemParam.LS, systemParam.A);
    if(check_cvode_flag(&flag, "CVodeSetLinearSolver", 1)) exit(-1);

    if(lsolver == 0) {
      /* Set the user-supplied Jacobian routine Jac */
      flag = CVodeSetJacFn(systemParam.cvodeMemPtr, jac_full_dense);
      if(check_cvode_flag(&flag, "CVDlsSetJacFn", 1)) exit(-1);
    } else if (lsolver == 1) {
      flag = CVodeSetJacFn(systemParam.cvodeMemPtr, NULL);
      if(check_cvode_flag(&flag, "CVDlsSetJacFn", 1)) exit(-1);
    }
#endif
  } else {
#if defined SUNDIALS2
    flag = CVSpgmr(systemParam.cvodeMemPtr, PREC_LEFT, 0);
    if(check_cvode_flag(&flag, "CVSpgmr", 1)) exit(-1);

    flag = CVSpilsSetGSType(systemParam.cvodeMemPtr, MODIFIED_GS);
    if(check_cvode_flag(&flag, "CVSpilsSetGSType", 1)) exit(-1);

    flag = CVSpilsSetPreconditioner(systemParam.cvodeMemPtr, jac_full_prec_setup,
                                    jac_full_prec_solveV3);
    if(check_cvode_flag(&flag, "CVSpilsSetPreconditioner", 1)) exit(-1);
#elif defined SUNDIALS3
    systemParam.LS = SUNSPGMR(y0, PREC_LEFT, 0);
    flag = CVSpilsSetLinearSolver(systemParam.cvodeMemPtr, systemParam.LS);
    if(check_cvode_flag(&flag, "CVSpilsSetLinearSolver", 1)) exit(-1);

    flag = SUNSPGMRSetGSType(systemParam.LS, MODIFIED_GS);
    if(check_cvode_flag(&flag, "SUNSPGMRSetGSType", 1)) exit(-1);

    flag = CVSpilsSetPreconditioner(systemParam.cvodeMemPtr, jac_full_prec_setup,
    				  jac_full_prec_solveV3);
    if(check_cvode_flag(&flag, "CVSpilsSetPreconditioner", 1)) exit(-1);
#elif defined SUNDIALS4
    systemParam.NLS = SUNNonlinSol_Newton(y0);
    flag = CVodeSetNonlinearSolver(systemParam.cvodeMemPtr, systemParam.NLS);
    if(check_cvode_flag(&flag, "CVodeSetNonlinearSolver", 1)) exit(-1);

    systemParam.LS = SUNLinSol_SPGMR(y0, PREC_LEFT, 0);
    flag = CVodeSetLinearSolver(systemParam.cvodeMemPtr, systemParam.LS, NULL);
    if(check_cvode_flag(&flag, "CVodeSetLinearSolver", 1)) exit(-1);

    flag = CVodeSetPreconditioner(systemParam.cvodeMemPtr, jac_full_prec_setup,
                                  jac_full_prec_solveV3);
    if(check_cvode_flag(&flag, "CVodeSetPreconditioner", 1)) exit(-1);
#endif
  }

  /* Set the maximum number of internal steps per CVode call and the maximum
   * allowable internal steps. */
  flag = CVodeSetMaxNumSteps(systemParam.cvodeMemPtr, maxsteps);
  if (check_cvode_flag(&flag, "CVodeSetMaxNumSteps", 1)) exit(-1);
  flag = CVodeSetMaxStep(systemParam.cvodeMemPtr, maxdt_internal);
  if (check_cvode_flag(&flag, "CVodeSetMaxStep", 1)) exit(-1);

  flag = CVodeSetNonlinConvCoef(systemParam.cvodeMemPtr, nlconvcoef);
  if (check_cvode_flag(&flag, "CVodeSetNonlinConvCoef", 1)) exit(-1);
#if defined SUNDIALS2 || defined SUNDIALS3
  flag = CVSpilsSetEpsLin(systemParam.cvodeMemPtr, epslin);
#else
  flag = CVodeSetEpsLin(systemParam.cvodeMemPtr, epslin);
#endif
  if (check_cvode_flag(&flag, "CVodeSetEpsLin", 1)) exit(-1);

  N_VDestroy_Serial(y0);

#ifdef USE_MPI
  systemParam.rank_weights =(double *)malloc(sizeof(double)*systemParam.nranks);
  double my_weight = 1.0;
  MPI_Allgather(&my_weight,1,MPI_DOUBLE,systemParam.rank_weights,1,MPI_DOUBLE,MPI_COMM_WORLD);
  double sum_weights = 0;
  for(int i = 0 ; i < systemParam.nranks; ++i)
  {
    sum_weights += systemParam.rank_weights[i];
  }
  double weight_factor = systemParam.nranks/sum_weights;
  for(int i = 0 ; i < systemParam.nranks; ++i)
  {
    systemParam.rank_weights[i] *= weight_factor;
  }
#endif
}

extern "C"
void zerork_cfd_plugin_setup_defaults
(
    double rtol,
    double atol,
    double sparsethresh
)
{
    int abstol_dens = 0; //mass fraction based abstol
    int maxsteps = 1e6;
    int maxord = 5;  //CVODE Default
    double maxdt_internal = DBL_MAX; //TODO: Reasonable value?
    double nlconvcoef = 0.1; //CVODE Default
    double epslin = 0.05;    //CVODE Default
    int multireac = 0;  //No GPU
    int gpu_id = -1;    //No GPU
    int lsolver = 2; //Sparse Iterative Preconditioned

    const char * ck_log_filename = "";       //No output
    const char * reactor_log_filename = "";  //No output

    int verbosity = 0;

    zerork_cfd_plugin_setup_full(verbosity,maxsteps,maxord,rtol,atol,abstol_dens,
                           maxdt_internal,sparsethresh,nlconvcoef,epslin,
                           lsolver,"mech.dat","therm.dat",ck_log_filename,
                           reactor_log_filename,&multireac,&gpu_id);
}


typedef std::pair<int,int> indexed_int;
typedef std::pair<int,double> indexed_double;
template<typename T, typename U>
struct lesser_second {
  bool operator()(std::pair<T,U> const &left,
                  std::pair<T,U> const &right)
  {
    return left.second < right.second;
  }
};
template<typename T, typename U>
struct greater_second {
  bool operator()(std::pair<T,U> const &left,
                  std::pair<T,U> const &right)
  {
    return left.second > right.second;
  }
};

extern "C"
void zerork_solve_reactors(const int nReactors, const double dt,
                        const double *T, const double *P,
                        double *massFracPtr,
                        double* reactorCost,
                        double* reactorGpu)
{
  static int n_calls = 0;
  n_calls++;

  //Timing info
  double startTime,stopTime,startAllTime,allTime,commTime,sortTime;
  commTime=0.0;
  sortTime=0.0;

  int nreactors_calc = nReactors;

  std::vector<int> reactors_on_rank(systemParam.nranks);
  reactors_on_rank[systemParam.rank] = nreactors_calc;

#ifdef USE_MPI

  MPI_Allgather(&nreactors_calc,1,MPI_INT,&reactors_on_rank[0],1,MPI_INT,MPI_COMM_WORLD);

#ifdef ZERORK_WORK_SHARING
  int max_recv = 0;
  int max_send = 0;

  //Counts of how many reactors we want to add to ours to be balanced overall
  std::vector<int> reactors_wanted(systemParam.nranks);
  std::vector<indexed_int> reactor_deficit(systemParam.nranks);

  //Sparse matrix of reactors transmitted to/from ranks.
  //  (N.B. Column indicies are in reverse order from normal CSR)
  std::vector<int> comm_mtx_row_sum(systemParam.nranks+1,0);
  std::vector<int> comm_mtx_col_idx;
  std::vector<int> comm_mtx_reactors;
  comm_mtx_col_idx.reserve(2*systemParam.nranks); //conservative estimate of final size
  comm_mtx_reactors.reserve(2*systemParam.nranks); //conservative estimate of final size

  int total_reactors = 0;
  int sorted_rank = -1;
  std::vector<int> rank_has_gpu(systemParam.nranks,0);
  int have_gpu = 0;
  double n_gpus = 0;
  if(systemParam.nranks>1)
  {
    startTime = getHighResolutionTime();

    //For weighting method 1.
    MPI_Allgather(&have_gpu,1,MPI_INT,&rank_has_gpu[0],1,MPI_INT,MPI_COMM_WORLD);

    double sum_weights = 0.0;
    for(int i = 0; i < systemParam.nranks; ++i)
    {
      total_reactors += reactors_on_rank[i];
      sum_weights += systemParam.rank_weights[i];
      n_gpus += rank_has_gpu[i];
    }

    int total_wanted = 0;
    double weight_factor = total_reactors/sum_weights;
    for(int i = 0; i < systemParam.nranks; ++i)
    {
      reactors_wanted[i] = systemParam.rank_weights[i]*weight_factor;
      total_wanted += reactors_wanted[i];
    }

    // Could have between 0 and nranks remaining reactors
    int reactor_remainders = total_reactors - total_wanted;
    assert(abs(reactor_remainders) < systemParam.nranks);
    int reactors_accounted = 0;
    for(int i = 0; i < systemParam.nranks; ++i)
    {
      if(i < abs(reactor_remainders))
      {
        if(reactor_remainders > 0) reactors_wanted[i] += 1;
        if(reactor_remainders < 0) reactors_wanted[i] -= 1;
      }
      reactor_deficit[i] = std::make_pair(i,reactors_wanted[i] - reactors_on_rank[i]);
      reactors_accounted += reactors_wanted[i];
    }
    assert(reactors_accounted == total_reactors);
    std::sort(reactor_deficit.begin(),reactor_deficit.end(),greater_second<int,int>());
    for(int i = 0; i < systemParam.nranks; ++i) { //note where we are in the list
      if(reactor_deficit[i].first == systemParam.rank) {
        sorted_rank = i;
        break;
      }
    }
    assert(sorted_rank >= 0);

    const int MIN_CHUNK = std::max(total_reactors/(systemParam.nranks*24),8);

    int give_rank = systemParam.nranks-1; // smallest deficit (i.e. negative deficits)
                              // gets stolen from first
    for(int i = 0; i < systemParam.nranks; ++i)
    {
      comm_mtx_row_sum[i+1] = comm_mtx_row_sum[i]; //keep the running sum going
      while(reactor_deficit[i].second >= MIN_CHUNK &&// while we stil want more
            give_rank > i)                        //and there's more to be taken
      {
        int extra_reactors_wanted    = reactor_deficit[i].second;
        int extra_reactors_available = -reactor_deficit[give_rank].second;
        if(extra_reactors_available >= MIN_CHUNK)
        {
          int reactors_taken = std::min(extra_reactors_available, extra_reactors_wanted);
          reactor_deficit[i].second -= reactors_taken;
          reactor_deficit[give_rank].second += reactors_taken;
          comm_mtx_row_sum[i+1] += 1;
          comm_mtx_col_idx.push_back(give_rank);
          comm_mtx_reactors.push_back(reactors_taken);
          //keep track of max storage requirement
          if(reactor_deficit[give_rank].first == systemParam.rank)
          { max_send = std::max(reactors_taken, max_send); }
          if(reactor_deficit[i].first == systemParam.rank)
          { max_recv = std::max(reactors_taken, max_recv); }
        }
        if(-reactor_deficit[give_rank].second < MIN_CHUNK)
        { //no more to give on this rank
          give_rank -= 1;
        }
        if(reactor_deficit[give_rank].second > 0)
        { //no more to give at all
          break;
        }
      }
    }
    nreactors_calc = reactors_wanted[systemParam.rank] - reactor_deficit[sorted_rank].second;
    int total_reactors_calc=0;
    MPI_Allreduce(&nreactors_calc,&total_reactors_calc,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    assert(total_reactors_calc == total_reactors);
    commTime = getHighResolutionTime() - startTime;
  }
#endif
#endif

  //Allocate space for reactors
  int nreactors_alloc = max(nReactors,nreactors_calc);
  std::vector<double> T_loc(nreactors_alloc);
  std::vector<double> P_loc(nreactors_alloc);
  std::vector<double> mf_loc(nreactors_alloc*systemParam.nSpc);
  std::vector<double> nstep_loc(nreactors_alloc);
  std::vector<double> gpu_flag_loc(nreactors_alloc);
  std::vector<double> cost_loc(nreactors_alloc);

  std::vector<double> reactorTime(nreactors_calc,0.0);
  std::vector<int> gpuSolve(nreactors_calc,0);

  bool sort_reactors = false;
  std::vector<indexed_double> dtEst;
  if(sort_reactors)
  {
    startTime = getHighResolutionTime();
    dtEst.resize(nreactors_alloc);
    for(int k = 0; k < nReactors; ++k)
    {
      if(reactorCost[k] == 0.0) reactorCost[k] = 1.0;
      dtEst[k] = std::make_pair(k,reactorCost[k]); //1/(number of steps to complete)
    }
    //sort(&dtEst[0],&dtEst[nReactors],lesser_second<int,double>()); //Smallest first
    sort(dtEst.begin(),dtEst.end(),greater_second<int,double>()); //Biggest first
    for(int k = 0; k < nReactors; ++k)
    {
      T_loc[dtEst[k].first] = T[k];
      P_loc[dtEst[k].first] = P[k];
      memcpy(&mf_loc[dtEst[k].first*systemParam.nSpc],
             &massFracPtr[k*systemParam.nSpc],
             sizeof(double)*systemParam.nSpc);
      nstep_loc[dtEst[k].first] = reactorCost[k];
      gpu_flag_loc[dtEst[k].first] = reactorGpu[k];
      cost_loc[dtEst[k].first] = reactorCost[k];
    }
    sortTime = getHighResolutionTime() - startTime;
  }
  else
  {
    //TODO: could remove this by using pointers instead of
    //      or in addition to std::vectors
    memcpy(&T_loc[0],T,sizeof(double)*nReactors);
    memcpy(&P_loc[0],P,sizeof(double)*nReactors);
    memcpy(&mf_loc[0],massFracPtr,sizeof(double)*nReactors*systemParam.nSpc);
    memcpy(&nstep_loc[0],reactorCost,sizeof(double)*nReactors);
    memcpy(&gpu_flag_loc[0],reactorGpu,sizeof(double)*nReactors);
    memcpy(&cost_loc[0],reactorCost,sizeof(double)*nReactors);
  }

#ifdef USE_MPI
#ifdef ZERORK_WORK_SHARING
  //Exchange excess reactors
  if(systemParam.nranks>1)
  {
    startTime = getHighResolutionTime();
    const int EXCHANGE_SEND_TAG = 42;
    if(nReactors < nreactors_calc)
    {
      MPI_Status status;
      int recv_idx = nReactors; //memory location to take in new reactors
      int recv_count_per_reactor = 5+systemParam.nSpc;
      if(sort_reactors) recv_count_per_reactor += 1;
      std::vector<double> recv_buf(recv_count_per_reactor*max_recv);
      for(int i = comm_mtx_row_sum[sorted_rank]; i < comm_mtx_row_sum[sorted_rank+1]; ++i)
      {
        //Bring 'em in
        int sorted_send_rank = comm_mtx_col_idx[i];
        int send_rank = reactor_deficit[sorted_send_rank].first;
        int recv_nreactors = comm_mtx_reactors[i];
        MPI_Recv(&recv_buf[0],recv_count_per_reactor*recv_nreactors,
                 MPI_DOUBLE, send_rank, EXCHANGE_SEND_TAG,MPI_COMM_WORLD,
                 &status);
        //Unpack 'em
        int buf_idx = 0;
        memcpy(&mf_loc[recv_idx*systemParam.nSpc], &recv_buf[buf_idx],
               sizeof(double)*systemParam.nSpc*recv_nreactors);
        buf_idx += recv_nreactors*systemParam.nSpc;
        memcpy(&T_loc[recv_idx], &recv_buf[buf_idx],
               sizeof(double)*recv_nreactors);
        buf_idx += recv_nreactors;
        memcpy(&P_loc[recv_idx], &recv_buf[buf_idx],
               sizeof(double)*recv_nreactors);
        buf_idx += recv_nreactors;
        memcpy(&nstep_loc[recv_idx], &recv_buf[buf_idx],
               sizeof(double)*recv_nreactors);
        buf_idx += recv_nreactors;
        memcpy(&gpu_flag_loc[recv_idx], &recv_buf[buf_idx],
               sizeof(double)*recv_nreactors);
        buf_idx += recv_nreactors;
        memcpy(&cost_loc[recv_idx], &recv_buf[buf_idx],
               sizeof(double)*recv_nreactors);
        if(sort_reactors) {
            buf_idx += recv_nreactors;
            for(int k = 0; k < recv_nreactors; ++k) {
                dtEst[recv_idx+k] = std::make_pair(recv_idx+k,recv_buf[buf_idx+k]);
            }
        }
        recv_idx += recv_nreactors;
      }
    }

    if(nReactors > nreactors_calc)
    {
      int send_idx = nReactors;
      int send_count_per_reactor = 5+systemParam.nSpc;
      if(sort_reactors) send_count_per_reactor += 1;
      std::vector<double> send_buf(max_send*send_count_per_reactor);
      for(int ircv = 0; ircv < systemParam.nranks; ++ircv)
      {
        int recv_rank = reactor_deficit[ircv].first;
        if(recv_rank == systemParam.rank) continue;
        for(int j = comm_mtx_row_sum[ircv]; j < comm_mtx_row_sum[ircv+1]; ++j)
        {
          int sorted_send_rank = comm_mtx_col_idx[j];
          int send_rank = reactor_deficit[sorted_send_rank].first;
          if(send_rank == systemParam.rank)
          {
            //Pack 'em up
            int send_nreactors = comm_mtx_reactors[j];
            send_idx -= send_nreactors;
            int buf_idx = 0;
            memcpy(&send_buf[buf_idx], &mf_loc[send_idx*systemParam.nSpc],
                   sizeof(double)*systemParam.nSpc*send_nreactors);
            buf_idx += send_nreactors*systemParam.nSpc;
            memcpy(&send_buf[buf_idx], &T_loc[send_idx],
                   sizeof(double)*send_nreactors);
            buf_idx += send_nreactors;
            memcpy(&send_buf[buf_idx], &P_loc[send_idx],
                   sizeof(double)*send_nreactors);
            buf_idx += send_nreactors;
            memcpy(&send_buf[buf_idx], &nstep_loc[send_idx],
                   sizeof(double)*send_nreactors);
            buf_idx += send_nreactors;
            memcpy(&send_buf[buf_idx], &gpu_flag_loc[send_idx],
                   sizeof(double)*send_nreactors);
            buf_idx += send_nreactors;
            memcpy(&send_buf[buf_idx], &cost_loc[send_idx],
                   sizeof(double)*send_nreactors);
            if(sort_reactors) {
                buf_idx += send_nreactors;
                for(int k = 0; k < send_nreactors; ++k) {
                    send_buf[buf_idx+k] = dtEst[send_idx+k].second;
                }
            }
            //Ship 'em out
            MPI_Send(&send_buf[0],send_count_per_reactor*send_nreactors,
                     MPI_DOUBLE, recv_rank, EXCHANGE_SEND_TAG, MPI_COMM_WORLD);
          }
        }
      }
    }
    commTime += getHighResolutionTime() - startTime;
  }

  //To maintain sorting we need to sort again
  if(sort_reactors && nreactors_calc > nReactors) {
    startTime = getHighResolutionTime();
    //sort(&dtEst[0],&dtEst[nreactors_calc],lesser_second<int,double>()); //Smallest first
    sort(dtEst.begin(),dtEst.end(),greater_second<int,double>()); //Biggest first
    std::vector<double> T_loc_tmp(nreactors_calc);
    std::vector<double> P_loc_tmp(nreactors_calc);
    std::vector<double> nstep_loc_tmp(nreactors_calc);
    std::vector<double> gpu_flag_loc_tmp(nreactors_calc);
    std::vector<double> cost_loc_tmp(nreactors_calc);
    std::vector<double> mf_loc_tmp(nreactors_calc*systemParam.nSpc);
    for(int k = 0; k < nReactors; ++k)
    {
      T_loc_tmp[dtEst[k].first] = T[k];
      P_loc_tmp[dtEst[k].first] = P[k];
      nstep_loc_tmp[dtEst[k].first] = reactorCost[k];
      gpu_flag_loc_tmp[dtEst[k].first] = reactorGpu[k];
      cost_loc_tmp[dtEst[k].first] = reactorCost[k];
      memcpy(&mf_loc_tmp[dtEst[k].first*systemParam.nSpc],
             &massFracPtr[k*systemParam.nSpc],
             sizeof(double)*systemParam.nSpc);
    }
    for(int k = nReactors; k < nreactors_calc; ++k)
    {
      T_loc_tmp[dtEst[k].first] = T_loc[k];
      P_loc_tmp[dtEst[k].first] = P_loc[k];
      nstep_loc_tmp[dtEst[k].first] = nstep_loc[k];
      gpu_flag_loc_tmp[dtEst[k].first] = gpu_flag_loc[k];
      cost_loc_tmp[dtEst[k].first] = cost_loc[k];
      memcpy(&mf_loc_tmp[dtEst[k].first*systemParam.nSpc],
             &mf_loc[k*systemParam.nSpc],
             sizeof(double)*systemParam.nSpc);
    }
    memcpy(&T_loc[0],&T_loc_tmp[0],sizeof(double)*nreactors_calc);
    memcpy(&P_loc[0],&P_loc_tmp[0],sizeof(double)*nreactors_calc);
    memcpy(&nstep_loc[0],&nstep_loc_tmp[0],sizeof(double)*nreactors_calc);
    memcpy(&gpu_flag_loc[0],&gpu_flag_loc_tmp[0],sizeof(double)*nreactors_calc);
    memcpy(&cost_loc[0],&cost_loc_tmp[0],sizeof(double)*nreactors_calc);
    memcpy(&mf_loc[0],&mf_loc_tmp[0],sizeof(double)*nreactors_calc*systemParam.nSpc);
    sortTime += getHighResolutionTime() - startTime;
  }
#endif
#endif

  startAllTime = getHighResolutionTime();

  //start with one step to avoid dbz in rare case that all reactors failed.
  double nstep_avg_gpu = 1.0;
  double nstep_avg_cpu = 1.0;
  int nGpuGroups = 0;
  {
    for(int k = 0; k < nreactors_calc; ++k)
    {
        startTime=getHighResolutionTime();
        nstep_loc[k] = zerork_cfd_plugin_solve(dt, T_loc[k], P_loc[k], &(mf_loc[k*systemParam.nSpc]));
        gpu_flag_loc[k] = 0.0;
        cost_loc[k] = nstep_loc[k];
        nstep_avg_cpu += nstep_loc[k];
        reactorTime[k] = getHighResolutionTime()-startTime;
    }
  }
  allTime = getHighResolutionTime() - startAllTime;

  //Unsort before sending back
  if(sort_reactors && nreactors_calc > nReactors)
  {
    startTime = getHighResolutionTime();
    std::vector<double> nstep_loc_tmp(nreactors_calc);
    std::vector<double> gpu_flag_loc_tmp(nreactors_calc);
    std::vector<double> cost_loc_tmp(nreactors_calc);
    std::vector<double> mf_loc_tmp(nreactors_calc*systemParam.nSpc);
    for(int k = 0; k < nreactors_calc; ++k)
    {
      nstep_loc_tmp[k] = nstep_loc[dtEst[k].first];
      gpu_flag_loc_tmp[k] = gpu_flag_loc[dtEst[k].first];
      cost_loc_tmp[k] = cost_loc[dtEst[k].first];
      memcpy(&mf_loc_tmp[k*systemParam.nSpc],
             &mf_loc[dtEst[k].first*systemParam.nSpc],
             sizeof(double)*systemParam.nSpc);
    }
    memcpy(&nstep_loc[0],&nstep_loc_tmp[0],sizeof(double)*nreactors_calc);
    memcpy(&gpu_flag_loc[0],&gpu_flag_loc_tmp[0],sizeof(double)*nreactors_calc);
    memcpy(&cost_loc[0],&cost_loc_tmp[0],sizeof(double)*nreactors_calc);
    memcpy(&mf_loc[0],&mf_loc_tmp[0],sizeof(double)*nreactors_calc*systemParam.nSpc);
    sortTime += getHighResolutionTime() - startTime;
  }

  //Per processor step counts and timing info
  int nGpuSolve = 0;
  int nCpuSolve = 0;
  double sumCpuReactorTime = 0.0;
  double sumGpuReactorTime = 0.0;
  for(int k = 0; k < nreactors_calc; ++k)
  {
    if(gpuSolve[k] > 0) {
        sumGpuReactorTime += reactorTime[k];
        ++nGpuSolve;
    } else {
        sumCpuReactorTime += reactorTime[k];
        ++nCpuSolve;
    }
  }
  double nstep_avg = nreactors_calc > 0 ? (nstep_avg_cpu + nstep_avg_gpu)/nreactors_calc: 0;
  nstep_avg_gpu = nGpuSolve > 0 ? (nstep_avg_gpu)/nGpuSolve : 0;
  nstep_avg_cpu = nCpuSolve > 0 ? (nstep_avg_cpu)/nCpuSolve : 0;

  int nReactorsTotal = nreactors_calc;
  double cpuPerStepTime = nCpuSolve > 0 ? sumCpuReactorTime/(nstep_avg_cpu*nCpuSolve) : 0;
  double gpuPerStepTime = nGpuSolve > 0 ? sumGpuReactorTime/(nstep_avg_gpu*nGpuSolve) : 0;

#ifdef USE_MPI
  startTime = getHighResolutionTime();
  MPI_Barrier(MPI_COMM_WORLD);
  double synchTime = getHighResolutionTime() - startTime;
  if(systemParam.nranks>1)
  {
#ifdef ZERORK_WORK_SHARING
    // Send back results
    startTime = getHighResolutionTime();
    const int EXCHANGE_RETURN_TAG = 43;
    if(nReactors < nreactors_calc)
    {
      int send_idx = nReactors;
      for(int i = comm_mtx_row_sum[sorted_rank]; i < comm_mtx_row_sum[sorted_rank+1]; ++i)
      {
        int sorted_recv_rank = comm_mtx_col_idx[i];
        int recv_rank = reactor_deficit[sorted_recv_rank].first;
        int send_nreactors = comm_mtx_reactors[i];
        MPI_Send(&nstep_loc[send_idx],
                 send_nreactors,
                 MPI_DOUBLE, recv_rank, EXCHANGE_RETURN_TAG,MPI_COMM_WORLD);
        MPI_Send(&gpu_flag_loc[send_idx],
                 send_nreactors,
                 MPI_DOUBLE, recv_rank, EXCHANGE_RETURN_TAG,MPI_COMM_WORLD);
        MPI_Send(&cost_loc[send_idx],
                 send_nreactors,
                 MPI_DOUBLE, recv_rank, EXCHANGE_RETURN_TAG,MPI_COMM_WORLD);
        MPI_Send(&mf_loc[send_idx*systemParam.nSpc],
                 systemParam.nSpc*send_nreactors,
                 MPI_DOUBLE, recv_rank, EXCHANGE_RETURN_TAG,MPI_COMM_WORLD);
        send_idx += send_nreactors;
      }
    }

    if(nReactors > nreactors_calc)
    {
      MPI_Status status;
      int recv_idx = nReactors;
      for(int ircv = 0; ircv < systemParam.nranks; ++ircv)
      {
        int send_rank = reactor_deficit[ircv].first;
        if(send_rank == systemParam.rank) continue;
        for(int j = comm_mtx_row_sum[ircv]; j < comm_mtx_row_sum[ircv+1]; ++j)
        {
          int sorted_recv_rank = comm_mtx_col_idx[j];
          int recv_rank = reactor_deficit[sorted_recv_rank].first;
          if(recv_rank == systemParam.rank)
          {
            int recv_nreactors = comm_mtx_reactors[j];
            recv_idx -= recv_nreactors;
            MPI_Recv(&nstep_loc[recv_idx],
                     recv_nreactors,
                     MPI_DOUBLE, send_rank, EXCHANGE_RETURN_TAG, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&gpu_flag_loc[recv_idx],
                     recv_nreactors,
                     MPI_DOUBLE, send_rank, EXCHANGE_RETURN_TAG, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&cost_loc[recv_idx],
                     recv_nreactors,
                     MPI_DOUBLE, send_rank, EXCHANGE_RETURN_TAG, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&mf_loc[recv_idx*systemParam.nSpc],
                     systemParam.nSpc*recv_nreactors,
                     MPI_DOUBLE, send_rank, EXCHANGE_RETURN_TAG, MPI_COMM_WORLD,
                     &status);
          }
        }
      }
    }
    commTime += getHighResolutionTime() - startTime;
#endif

    startTime = getHighResolutionTime();
    //Output info on reactors calculated per rank
    std::vector<int> nreactors_calc_ranks(systemParam.nranks);
    std::vector<double> all_time_ranks(systemParam.nranks);
    MPI_Gather(&nreactors_calc,1,MPI_INT,&nreactors_calc_ranks[0],1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&allTime,1,MPI_DOUBLE,&all_time_ranks[0],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(systemParam.verbosity > 0) {
      if(systemParam.rank==0) {
        double total_time = 0.0;
        double max_time = 0.0;
        for(int i = 0; i < systemParam.nranks; ++i) {
          printf("Rank %d calculated %d reactors in %f seconds.\n",i,nreactors_calc_ranks[i],all_time_ranks[i]);
          total_time += all_time_ranks[i];
          max_time = std::max(max_time,all_time_ranks[i]);
        }
        double avg_time = total_time/systemParam.nranks;
        double wasted_time = max_time - avg_time;
        printf("Max Time, Avg Time, Wasted Time = %f, %f, %f\n",max_time,avg_time,wasted_time);
      }
    }
    commTime += getHighResolutionTime() - startTime;

#ifdef ZERORK_WORK_SHARING
    startTime = getHighResolutionTime();
    //Update rank_weights
    if(systemParam.rank==0 && systemParam.nranks > 1) {
      const int update_method = 1;
      if(update_method == -1)
      {
        //no update
      }
      else if(update_method == 0)
      {
        if(n_gpus == 0 || n_gpus == systemParam.nranks)
        {
          for(int i = 0; i < systemParam.nranks; ++i)
          {
            systemParam.rank_weights[i] = 1.0;
          }
        }
        else
        {
          double total_all_time_ranks = 0.0;
          double total_all_time_gpus = 0.0;
          double total_all_time_cpus = 0.0;

          double prev_gpu_weight = 0.0;
          double prev_cpu_weight = 0.0;

          int nreactors_calc_gpu = 0;
          int nreactors_calc_cpu = 0;

          for(int i = 0; i < systemParam.nranks; ++i)
          {
            total_all_time_ranks += all_time_ranks[i];
            if(rank_has_gpu[i] == 1)
            {
              total_all_time_gpus += all_time_ranks[i];
              prev_gpu_weight = systemParam.rank_weights[i];
              nreactors_calc_gpu = nreactors_calc_ranks[i];
            }
            else
            {
              total_all_time_cpus += all_time_ranks[i];
              prev_cpu_weight = systemParam.rank_weights[i];
              nreactors_calc_cpu = nreactors_calc_ranks[i];
            }
          }
          const double alpha = 0.8; //Don't over correct
          double sum_weights = 0.0;
          double new_gpu_weight=prev_gpu_weight;
          double new_cpu_weight=prev_cpu_weight;

          if(total_all_time_gpus > 0.0)
          {
              new_gpu_weight = ((double)nreactors_calc_gpu/(double)total_reactors)*
                               (total_all_time_ranks/total_all_time_gpus)*
                               ((double)n_gpus/(double)systemParam.nranks);
          }
          if(total_all_time_cpus > 0.0)
          {
              new_cpu_weight = ((double)nreactors_calc_cpu/(double)total_reactors)*
                               (total_all_time_ranks/total_all_time_cpus)*
                               ((double)(systemParam.nranks-n_gpus)/(double)systemParam.nranks);
          }
          new_gpu_weight = (1.0-alpha)*prev_gpu_weight + alpha*new_gpu_weight;
          new_cpu_weight = (1.0-alpha)*prev_cpu_weight + alpha*new_cpu_weight;
          for(int i = 0; i < systemParam.nranks; ++i) {
            if(rank_has_gpu[i])
            {
               systemParam.rank_weights[i] = new_gpu_weight;
            }
            else
            {
               systemParam.rank_weights[i] = new_cpu_weight;
            }
            sum_weights += systemParam.rank_weights[i];
          }
          double weight_factor = systemParam.nranks/sum_weights;
          for(int i = 0 ; i < systemParam.nranks; ++i)
          {
            systemParam.rank_weights[i] *= weight_factor;
            if(systemParam.verbosity >= 2) {
              printf("RANK[%d] weight: %f\n",i,systemParam.rank_weights[i]);
            }
          }
        }
      }
      else if(update_method == 1)
      {
        double total_all_time_ranks = 0.0;
        for(int i = 0; i < systemParam.nranks; ++i) {
          total_all_time_ranks += all_time_ranks[i];
        }
        const double alpha = 0.1; //Don't over correct
        double sum_weights = 0.0;
        for(int i = 0; i < systemParam.nranks; ++i) {
          double last_weight = systemParam.rank_weights[i];
          double new_weight;
          if(nreactors_calc_ranks[i] > 0)
          {
              new_weight = ((double)nreactors_calc_ranks[i]/(double)total_reactors)*
                                        (total_all_time_ranks/all_time_ranks[i]);
          }
          else
          {
            new_weight = last_weight;
          }
          systemParam.rank_weights[i] =(1.0-alpha)*last_weight + alpha*new_weight;
          sum_weights += systemParam.rank_weights[i];
        }
        double weight_factor = systemParam.nranks/sum_weights;
        for(int i = 0 ; i < systemParam.nranks; ++i)
        {
          systemParam.rank_weights[i] *= weight_factor;
          if(systemParam.verbosity >= 2) {
            printf("RANK[%d] weight: %f\n",i,systemParam.rank_weights[i]);
          }
        }
      }
      else
      {
        printf("Invalid weighting update method %d.  Quitting.\n",update_method);
        exit(-1);
      }
    }
    MPI_Bcast(systemParam.rank_weights,systemParam.nranks,MPI_DOUBLE,0,MPI_COMM_WORLD);
    commTime += getHighResolutionTime() - startTime;
#endif
  }
#endif //USE_MPI

  //Unsort now that we have all ours back
  if(sort_reactors && nreactors_calc <= nReactors)
  {
    startTime = getHighResolutionTime();
    std::vector<double> nstep_loc_tmp(nReactors);
    std::vector<double> gpu_flag_loc_tmp(nReactors);
    std::vector<double> cost_loc_tmp(nReactors);
    std::vector<double> mf_loc_tmp(nReactors*systemParam.nSpc);
    for(int k = 0; k < nReactors; ++k)
    {
      nstep_loc_tmp[k] = nstep_loc[dtEst[k].first];
      gpu_flag_loc_tmp[k] = gpu_flag_loc[dtEst[k].first];
      cost_loc_tmp[k] = cost_loc[dtEst[k].first];
      memcpy(&mf_loc_tmp[k*systemParam.nSpc],
             &mf_loc[dtEst[k].first*systemParam.nSpc],
             sizeof(double)*systemParam.nSpc);
    }
    memcpy(&nstep_loc[0],&nstep_loc_tmp[0],sizeof(double)*nReactors);
    memcpy(&gpu_flag_loc[0],&gpu_flag_loc_tmp[0],sizeof(double)*nReactors);
    memcpy(&cost_loc[0],&cost_loc_tmp[0],sizeof(double)*nReactors);
    memcpy(&mf_loc[0],&mf_loc_tmp[0],sizeof(double)*nReactors*systemParam.nSpc);
    sortTime += getHighResolutionTime() - startTime;
  }

  double nstep_avg_global = nstep_avg;
#ifdef USE_MPI
  if(systemParam.nranks>1)
  {
    // Collect timing/step count data
    double rr; //reduced real
    int ri; //reduced int

    //Re-average after reduction
    nstep_avg_global *= nreactors_calc;
    nstep_avg_cpu *= nCpuSolve;
    nstep_avg_gpu *= nGpuSolve;

    MPI_Reduce(&n_calls,&ri,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) n_calls = ri;

    nReactorsTotal=nreactors_calc;
    MPI_Reduce(&nReactorsTotal,&ri,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) nReactorsTotal = ri;
    MPI_Reduce(&nCpuSolve,&ri,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) nCpuSolve= ri;
    MPI_Reduce(&nGpuSolve,&ri,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) nGpuSolve= ri;
    MPI_Reduce(&nGpuGroups,&ri,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) nGpuGroups = ri;

    MPI_Reduce(&nstep_avg_global,&rr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) nstep_avg_global = nReactorsTotal > 0 ? rr/nReactorsTotal : 0;
    MPI_Reduce(&nstep_avg_cpu,&rr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) nstep_avg_cpu = nCpuSolve > 0 ? rr/nCpuSolve : 0;
    MPI_Reduce(&nstep_avg_gpu,&rr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) nstep_avg_gpu = nGpuSolve > 0 ? rr/nGpuSolve : 0;

    //Calc per step times based on sum of times
    MPI_Reduce(&sumCpuReactorTime,&rr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) cpuPerStepTime = nCpuSolve > 0 ? rr/(nstep_avg_cpu*nCpuSolve) : 0;
    MPI_Reduce(&sumGpuReactorTime,&rr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) gpuPerStepTime = nGpuSolve > 0 ? rr/(nstep_avg_gpu*nGpuSolve) : 0;

    //Output times based on max times
    MPI_Reduce(&sumCpuReactorTime,&rr,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) sumCpuReactorTime = rr;
    MPI_Reduce(&sumGpuReactorTime,&rr,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) sumGpuReactorTime = rr;
    MPI_Reduce(&commTime,&rr,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) commTime = rr;
    MPI_Reduce(&synchTime,&rr,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) synchTime= rr;
    MPI_Reduce(&sortTime,&rr,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) sortTime= rr;
    MPI_Reduce(&allTime,&rr,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(systemParam.rank==0) allTime = rr;

    if(systemParam.rank==0 && systemParam.verbosity >= 2)
    {
      printf("commTime,synchTime: %f, %f\n",commTime, synchTime);
    }
  }
#endif
  if(systemParam.rank==0 && sort_reactors && systemParam.verbosity >= 2)
  {
    printf("sortTime:  %f\n", sortTime);
  }

  //Put our local mass fracs back into the array that was passed in.
  memcpy(massFracPtr,&mf_loc[0],sizeof(double)*systemParam.nSpc*nReactors);
  memcpy(reactorGpu,&gpu_flag_loc[0],sizeof(double)*nReactors);
  memcpy(reactorCost,&cost_loc[0],sizeof(double)*nReactors);

  double sumFuncTime = systemParam.colPermTime + systemParam.jacSetupTime +
                systemParam.precSetupTime + systemParam.jacFactorTime +
                systemParam.backsolveTime + systemParam.funcTime;
  double otherTime = sumCpuReactorTime + sumGpuReactorTime - sumFuncTime;

  if(systemParam.verbosity > 4)
  {
      //N.B. We don't have final temp or pressure just mass fracs
      for(int k = 0; k < nReactors; ++k)
      {
            printf( "%8d  %14.7e  %14.7e  %2d\n",
    	              k,massFracPtr[k*systemParam.nSpc],reactorTime[k],gpuSolve[k]);
      }
  }

  if(systemParam.verbosity > 0)
  {
//      fprintf(reactor_log_file,"nReactors, Min, Med, Max, Avg ="
//             "%d, %20.12g, %20.12g, %20.12g, %20.12g\n",
//             nReactors,
//             dtEst[0].second,
//             dtEst[(nReactors-1)/2].second,
//             dtEst[nReactors-1].second,
//             dt/nstep_avg);

//      printf("%15s%15s%15s%15s%15s%15s%15s%15s\n","cpTime","jsTime",
//             "psTime","jfTime","bsTime","fnTime","otTime","simTime");
//      printf("%14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g %14.10g\n",
//             systemParam.colPermTime,systemParam.jacSetupTime,
//             systemParam.precSetupTime,systemParam.jacFactorTime,
//             systemParam.backsolveTime,systemParam.funcTime,otherTime,sumFuncTime+otherTime);

//      fprintf(reactor_log_file,"\nreactor split - %d %d\n",nCpuSolve,nGpuSolve);
//      fprintf(reactor_log_file,"\ntotal time - %14.7e %14.7e %14.7e\n",sumCpuReactorTime,sumGpuReactorTime,allTime);
//      fprintf(reactor_log_file,"#%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s\n",
//                              "cells_solved","n_cpu","n_gpu","n_gpu_groups,"n_steps_avg",
//                              "n_steps_avg_cpu","n_steps_avg_gpu","time_cpu",
//                              "time_gpu","time_per_step_cpu","time_per_gpu_step",
//                              "time_total");

      if(systemParam.rank == 0)
      {
          fprintf(reactor_log_file," %12d",n_calls);
          fprintf(reactor_log_file," %16d %16d %16d %16d",
                                   nReactorsTotal,nCpuSolve,nGpuSolve,nGpuGroups);
          fprintf(reactor_log_file," %16.5e %16.5e %16.5e %16.5e",
                                   nstep_avg_global,nstep_avg_cpu,
                                   nstep_avg_gpu,sumCpuReactorTime);
          fprintf(reactor_log_file," %16.5e %16.5e %16.5e %16.5e\n",
                                  sumGpuReactorTime,cpuPerStepTime,
                                  gpuPerStepTime, allTime);
          fflush(reactor_log_file);
      }
  }

  if(systemParam.verbosity > 2) zerork_print_stats();
}


long int zerork_cfd_plugin_solve(const double dt, const double T,
                           const double P, double *massFracPtr)
{
  int flag;
  double tnext,tcurr;

  // reset the mass fraction, density and relative vol at new temperature
  systemParam.Press = P;
  systemParam.invDens = systemParam.mech->getDensityFromTPY(T,P,massFracPtr);
  systemParam.invDens=1.0/systemParam.invDens;

  double startTime = getHighResolutionTime();
  N_Vector systemState = N_VNew_Serial(systemParam.nSpc+1);
  memcpy(NV_DATA_S(systemState),massFracPtr,sizeof(double)*systemParam.nSpc);
  NV_Ith_S(systemState,systemParam.nSpc) = T;

  // reset the time
  tcurr=0.0;
  tnext=dt;

  // reinitialize cvode
  flag = CVodeReInit(systemParam.cvodeMemPtr, tcurr, systemState);

  // re-set tolerances
  if(systemParam.abstol_dens == 1) //use molar density based absolute tolerance
  {
      N_Vector atol_vector = N_VNew_Serial(systemParam.nSpc+1);
      double* atol_vector_ptr = NV_DATA_S(atol_vector);
      double reactor_density = 1.0/systemParam.invDens;
      for(int j=0;j<systemParam.nSpc;++j)
      {
          double molar_density =
              reactor_density*systemParam.invMolWt[j]*1.0e-3; //mks->cgs
          atol_vector_ptr[j] =
              systemParam.atol/molar_density;
      }
      atol_vector_ptr[systemParam.nSpc] = systemParam.atol; //Temperature stays the same
      flag = CVodeSVtolerances(systemParam.cvodeMemPtr,
                               systemParam.rtol, atol_vector);
      if (check_cvode_flag(&flag, "CVodeSVtolerances", 1)) exit(-1);
      N_VDestroy_Serial(atol_vector);
  }
  else
  {
      flag = CVodeSStolerances(systemParam.cvodeMemPtr,
                               systemParam.rtol, systemParam.atol);
      if (check_cvode_flag(&flag, "CVodeSStolerances", 1)) exit(-1);
  }

  flag = CVodeSetMaxStep(systemParam.cvodeMemPtr,
                         std::min(systemParam.maxdt_internal,tnext));

  int maxRetry = 5;
  int nRetry   = 0;
  while(1)
  {
      flag = CVode(systemParam.cvodeMemPtr, tnext, systemState, &tcurr, CV_NORMAL);

      if (check_cvode_flag(&flag, "CVode", 1)) {
        printf("WARNING: Failed integration current temperature: %g\n",NV_Ith_S(systemState,systemParam.nSpc));
        if(nRetry < maxRetry) {
          ++nRetry;
          CVodeReInit(systemParam.cvodeMemPtr, tcurr, systemState);
          continue;
        } else {
          break;
        }
      }
      if (flag == CV_SUCCESS) break;
  }

  zerork_update_cvode_stats(&systemParam, systemParam.cvodeMemPtr, false);
  systemParam.simTime += getHighResolutionTime() - startTime;
  long int nsteps;
  CVodeGetNumSteps(systemParam.cvodeMemPtr,&nsteps);

  //Failed to complete
  if(flag < 0) {
    printf("WARNING: Failed to complete CPU integration.\n");
    printf("         Check results carefully.\n");
  }

  memcpy(massFracPtr,NV_DATA_S(systemState),sizeof(double)*systemParam.nSpc);

  N_VDestroy_Serial(systemState);
  return nsteps;
}


extern "C" void zerork_cfd_plugin_close(void)
{
  // Free integrator memory
  CVodeFree(&(systemParam.cvodeMemPtr));

#if defined SUNDIALS3 || defined SUNDIALS4
  /* Free the linear solver memory */
  SUNLinSolFree(systemParam.LS);

  if(systemParam.lsolver < 2) {
    /* Free the matrix memory */
    SUNMatDestroy(systemParam.A);
  }
#endif
#if defined SUNDIALS4
  /* Free the linear solver memory */
  SUNNonlinSolFree(systemParam.NLS);
#endif

  // Free parameter memory
  free(systemParam.molWt);
  free(systemParam.invMolWt);
  free(systemParam.netProd);
  free(systemParam.fwdROP);
  free(systemParam.createRate);
  free(systemParam.destroyRate);
  free(systemParam.conc);
  free(systemParam.Energy);
  free(systemParam.CvMass);

  free_Jsparse(systemParam.sparseMtx);
  delete systemParam.mech;

#ifdef USE_MPI
  free(systemParam.rank_weights);
#endif

  // clean up procedures
  fclose(reactor_log_file);

  fflush(stdout);
}

static FILE* open_reactor_log_file(const char *filename)
{
  FILE * file_ptr=fopen(filename,"w");
  if(file_ptr==NULL) {
      printf("ERROR: could not open output file %s for write\n",
             filename);
      exit(-1);
  }
  return file_ptr;
}

void zerork_print_stats()
{
  static int print_count=0;
  if(print_count==0)
  {
    //print header
    printf("   nprnt");
    printf("   nstep   nertf  nnlsit  nnlcvf ngevals prsetup ncolprm nfactor"
           " backslv  intfcl     nnz    rnnz   lunnz    jac setup tm"
           "   prec setup tm     col perm tm      jac fac tm    backsolve tm"
           "     rhs func tm        other tm      sim tm\n");
  }

  double otherTime = systemParam.simTime -
                         ( systemParam.jacSetupTime + systemParam.precSetupTime
                          + systemParam.colPermTime + systemParam.jacFactorTime
                          + systemParam.backsolveTime + systemParam.funcTime);
  printf("  %6d",print_count);
  printf("  %6ld  %6ld  %6ld  %6ld  %6ld  %6d  %6d  %6d  %6d  %6d  %6d  %6d  %6d",
	     systemParam.nsteps,systemParam.netfails,
	     systemParam.nniters,systemParam.nncfails,
             systemParam.ngevals,systemParam.nJacSetup,
             systemParam.nColPerm,systemParam.nJacFactor,
             systemParam.nBackSolve,systemParam.nFunc,
	     systemParam.sparseMtx->nNonZero,
	     systemParam.sparseMtx->reduceNNZ,
	     systemParam.sparseMtx->LUnnz);
  printf("  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %8.4e\n",
             systemParam.jacSetupTime,
             systemParam.precSetupTime,
             systemParam.colPermTime,
             systemParam.jacFactorTime,
             systemParam.backsolveTime,
             systemParam.funcTime, otherTime, systemParam.simTime);
  ++print_count;
}

static void zerork_update_cvode_stats(cv_param* cvp, void* cvode_mem_ptr, bool multireac)
{
  int flag;
  long int nsteps,nfevals,nlinsetups,netfails;
  long int nniters,nncfails,ngevals;
  int qlast, qcur;
  long int nfeDQ = 0;
  double hinused,hlast,hcur,tcur;
  flag = CVodeGetIntegratorStats(cvode_mem_ptr, &nsteps, &nfevals,
                        &nlinsetups, &netfails, &qlast, &qcur,
                        &hinused, &hlast, &hcur, &tcur);
  flag = CVodeGetNonlinSolvStats(cvode_mem_ptr, &nniters, &nncfails);
  flag = CVodeGetNumGEvals(cvode_mem_ptr, &ngevals);
  cvp->nsteps += nsteps;
  cvp->nfevals += nfevals;
  cvp->nlinsetups += nlinsetups;
  cvp->netfails += netfails;
  cvp->nniters += nniters;
  cvp->nncfails += nncfails;
  cvp->ngevals += ngevals;
}

void zerork_reset_stats()
{
  systemParam.nFunc         = 0;
  systemParam.nJacSetup     = 0;
  systemParam.nJacFactor    = 0;
  systemParam.nBackSolve    = 0;
  systemParam.nJacRescale   = 0;
  systemParam.nColPerm      = 0;
  systemParam.colPermTime   = 0.0;
  systemParam.jacFactorTime = 0.0;
  systemParam.backsolveTime = 0.0;
  systemParam.jacSetupTime  = 0.0;
  systemParam.precSetupTime = 0.0;
  systemParam.funcTime      = 0.0;
  systemParam.simTime       = 0.0;
  systemParam.nsteps        = 0;
  systemParam.nfevals       = 0;
  systemParam.nlinsetups    = 0;
  systemParam.netfails      = 0;
  systemParam.nniters       = 0;
  systemParam.nncfails      = 0;
  systemParam.ngevals       = 0;
}

