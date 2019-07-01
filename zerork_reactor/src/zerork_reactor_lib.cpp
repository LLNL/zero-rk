
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

//CVIDT Code
#include "zerork_reactor_lib.h"
#include "cv_param_sparse.h"
#include "ode_funcs.h"
#include "matrix_funcs.h"
#include "utility_funcs.h"
#include "multizone.h"

//ZERORK

//Sundials
#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
#include <cvode/cvode_dense.h>      // prototypes & constants for CVDense
#include <cvode/cvode_spgmr.h>      // prototypes & constants for CVSPGMR

#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros
#include "ext/cvode_user_sparselu.h" // CVUserSuperLU
#include "ext/cvode_user_timed_dense.h" // CVUserTimedLapackDense


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
static double calc_max_dt_ratio(double dt);
static void zerork_update_cvode_stats(cv_param* cvp, void* cvode_mem_ptr, bool multireac);
int const_dpdt_wsr_sage(realtype t, N_Vector y, N_Vector ydot, void *user_data);

//Used by LLNL converge udf
extern "C"
void zerork_reactor_setup_full
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
    int CVCP,
    const char* mech_filename,
    const char* therm_filename,
    const char* ck_log_filename,
    const char* reactor_log_filename,
    int n_reactors_max,
    int n_reactors_min,
    double logiA, double logiK, double logiQ,
    double logiM, double logiB, double logirNu,
    int cusolver_factor_alg, int cusolver_solve_alg,
    int* multireac,
    int* gpu_id
)
{
  int j;
  int n_gpus = 0;
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

  systemParam.logiA = logiA;
  systemParam.logiK = logiK;
  systemParam.logiQ = logiQ;
  systemParam.logiM = logiM;
  systemParam.logiB = logiB;
  systemParam.logirNu = logirNu;

  *gpu_id = -1;
  *multireac = 0;

  int zerork_verbosity = 0;
  if(systemParam.rank==0) zerork_verbosity = 1;
  // Build zerork::mechanism
  if(*multireac)
  {
    std::cerr << " Multi-reactor disabled in non-cuda lib. QUITTING."
        << std::endl;
    exit(-1);
  }
  else
  {
      gasMech = new zerork::mechanism(mech_filename,therm_filename,
                                   ck_log_filename,zerork_verbosity);
  }

  nSpc=gasMech->getNumSpecies();
  nStep=gasMech->getNumSteps();
  nState=nSpc+1;

  if(*multireac)
  {
      systemParam.nReactors = n_reactors_max;
  }
  else
  {
      systemParam.nReactors = 1;
  }

  // set up the system parameters
  systemParam.maxsteps=maxsteps;
  systemParam.maxord=maxord;
  systemParam.rtol=rtol;
  systemParam.atol=atol;
  systemParam.abstol_dens = abstol_dens;
  systemParam.maxdt_internal=maxdt_internal;
  systemParam.nlconvcoef=nlconvcoef;
  systemParam.epslin=epslin;
  systemParam.gpu_id = *gpu_id;
  systemParam.verbosity = verbosity;
  systemParam.lsolver = lsolver;
  systemParam.nReactorsMax = n_reactors_max;
  systemParam.nReactorsMin = n_reactors_min;

  systemParam.doingJacSetup = false;
  systemParam.nSpc=nSpc;
  systemParam.minMassFrac=1.0e-30;
  systemParam.sqrtUnitRnd=sqrt(UNIT_ROUNDOFF);
  systemParam.Tref=1.0;  //TODO:Make functions use this value so it is adjustable
  systemParam.mech=gasMech;
  systemParam.Press=(double *)malloc(sizeof(double)*systemParam.nReactors);
  systemParam.invDens =(double *)malloc(sizeof(double)*systemParam.nReactors);
  systemParam.meanCvMass = (double *)malloc(sizeof(double)*systemParam.nReactors);
  systemParam.dpdt = (double *)malloc(sizeof(double)*systemParam.nReactors);
  systemParam.dTemp_dt = (double *)malloc(sizeof(double)*systemParam.nReactors);

  {
      systemParam.systemState_host=(double *)malloc(sizeof(double)*nState*systemParam.nReactors);
      systemParam.systemDeriv_host=(double *)malloc(sizeof(double)*nState*systemParam.nReactors);
      systemParam.tmp1_host=(double *)malloc(sizeof(double)*nState*systemParam.nReactors);
      systemParam.tmp2_host=(double *)malloc(sizeof(double)*nState*systemParam.nReactors);
      systemParam.tmp3_host=(double *)malloc(sizeof(double)*nState*systemParam.nReactors);
  }

  systemParam.Jac_host = NULL;

  systemParam.netProd    =(double *)malloc(sizeof(double)*nSpc*systemParam.nReactors);
  systemParam.Energy     =(double *)malloc(sizeof(double)*nSpc*systemParam.nReactors);
  systemParam.CvMass     =(double *)malloc(sizeof(double)*nSpc*systemParam.nReactors);
  systemParam.molWt      =(double *)malloc(sizeof(double)*nSpc);
  systemParam.invMolWt   =(double *)malloc(sizeof(double)*nSpc);
  systemParam.fwdROP     =(double *)malloc(sizeof(double)*nStep*systemParam.nReactors);
  systemParam.createRate =(double *)malloc(sizeof(double)*nSpc*systemParam.nReactors);
  systemParam.destroyRate=(double *)malloc(sizeof(double)*nSpc*systemParam.nReactors);
  systemParam.conc       =(double *)malloc(sizeof(double)*nSpc*systemParam.nReactors);

  zerork_reset_stats();


  // set constant parameters
  gasMech->getMolWtSpc(systemParam.molWt);
  for(j=0; j<nSpc; j++)
    {systemParam.invMolWt[j]=1.0/systemParam.molWt[j];}


  systemParam.sparseMtx=
        (Jsparse *)alloc_Jsparse(*(systemParam.mech),systemParam.nReactors);



  systemParam.nMatrixReactors = systemParam.nReactors;

  systemParam.nMatrixThreads = 1;
  const char* temp = getenv("ZERORK_NUM_MATRIX_THREADS");
  if(temp != NULL)
  {
     std::string env_str(temp);
     std::stringstream env_ss(env_str);
     env_ss >> systemParam.nMatrixThreads;
     if(!env_ss) //couldn't convert
     {
        std::cout << "Warning: Couldn't convert ZERORK_NUM_MATRIX_THREADS "
                  << "to integer.  Using value of 1." << std::endl;
        systemParam.nMatrixThreads = 1;
     }
     std::cout << "ZERORK_NUM_MATRIX_THREADS=" << systemParam.nMatrixThreads<< std::endl;
  }

  systemParam.sparseMtx->permThresh = 0.3;
  systemParam.sparseMtx->maxGammaChangeOrder = 3.0;
  systemParam.sparseMtx->strictSamePattern = 0;

  systemParam.sparseMtx->threshType = 1;
  systemParam.sparseMtx->offDiagThreshold = sparsethresh;
  systemParam.sparseMtx->ILU = false;
  systemParam.sparseMtx->fakeUpdate = false;

  for(j=0; j<systemParam.nReactors; ++j)
  {
    systemParam.sparseMtx->optionSLU[j].DiagPivotThresh = 0.0;
  }

  systemParam.sparseMtx->permutationType = 1;

  systemParam.sparseMtx->nReactors= systemParam.nReactors;
  for(j=0;j<systemParam.nReactors;++j)
  {
      systemParam.sparseMtx->reduceNNZ[j]=0; // won't be set until the first J
      systemParam.sparseMtx->LUnnz[j] = 0; // won't be set until the first J
      systemParam.sparseMtx->fillFactor[j]= 100.; // won't be set until the first J
      systemParam.sparseMtx->isFirstFactor[j] = 1;
      systemParam.sparseMtx->numPermReUses[j] = 0;
  }
#ifdef ZERORK_USE_CUSOLVERRF
  systemParam.sparseMtx->isFirstFactor_cusolverRf = 1;
#endif
#ifdef ZERORK_USE_CUSOLVERSP
  systemParam.sparseMtx->isFirstFactor_cusolverSp = 1;
#endif
  systemParam.prevNumErrTestFails = 0;


  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula and the use of a Newton iteration */
  systemParam.cvodeMemPtr = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_cvode_flag(systemParam.cvodeMemPtr, "CVodeCreate", 0)) exit(-1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */

  systemParam.constPress = CVCP == 1 ? true : false;

  N_Vector y0 = N_VNew_Serial(nState);
  if(systemParam.constPress)
  {
#ifdef ZERORK_USE_SAGE_RHS
      flag=CVodeInit(systemParam.cvodeMemPtr, const_dpdt_wsr_sage, 0.0, y0);
#else
      flag=CVodeInit(systemParam.cvodeMemPtr, const_dpdt_wsr, 0.0, y0);
#endif
  }
  else
  {
      flag=CVodeInit(systemParam.cvodeMemPtr, const_vol_wsr, 0.0, y0);
  }

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
    //flag = CVLapackDense(systemParam.cvodeMemPtr, nState);
    flag = CVUserTimedLapackDense(systemParam.cvodeMemPtr, nState);
    if(check_cvode_flag(&flag, "CVLapackDense", 1)) exit(-1);
    if(lsolver == 0) {
      flag = CVDlsSetDenseJacFn(systemParam.cvodeMemPtr, jac_full_dense);
      if(check_cvode_flag(&flag, "CVDlsDenseJacFn", 1)) exit(-1);
    } else if (lsolver == 1) {
      flag = CVDlsSetDenseJacFn(systemParam.cvodeMemPtr, NULL);
      if(check_cvode_flag(&flag, "CVDlsDenseJacFn", 1)) exit(-1);
    }
  } else if (lsolver == 2) {
    flag = CVUserSuperLU(systemParam.cvodeMemPtr,
                         nState,
                         systemParam.sparseMtx->nNonZero,
                         systemParam.sparseMtx->mtxRowIdx,
                         systemParam.sparseMtx->mtxColSum,
                         jac_full_sparse);
    if(check_cvode_flag(&flag, "CVUserSuperLU", 1)) exit(-1);
  } else if (lsolver == 3) {
    flag = CVUserSuperLU(systemParam.cvodeMemPtr,
                         nState,
                         systemParam.sparseMtx->nNonZero,
                         systemParam.sparseMtx->mtxRowIdx,
                         systemParam.sparseMtx->mtxColSum,
                         jac_full_sparse_divided_diff);
    if(check_cvode_flag(&flag, "CVUserSuperLU", 1)) exit(-1);
  } else { //default for lsolver != 0,1,2,3
    /* Call CVSpgmr to specify the linear solver CVSPGMR
       with left preconditioning and the maximum Krylov dimension maxl */
    flag = CVSpgmr(systemParam.cvodeMemPtr, PREC_LEFT, 0);
    if(check_cvode_flag(&flag, "CVSpgmr", 1)) exit(-1);

    /* Set modified Gram-Schmidt orthogonalization, preconditioner
       setup and solve routines Precond and PSolve, and the pointer
       to the user-defined block data */
    flag = CVSpilsSetGSType(systemParam.cvodeMemPtr, MODIFIED_GS);
    if(check_cvode_flag(&flag, "CVSpilsSetGSType", 1)) exit(-1);

    /* Set preconditioner setup and solve routines Precond and PSolve,
       and the pointer to the user-defined block data */
    flag = CVSpilsSetPreconditioner(systemParam.cvodeMemPtr,
                                    jac_full_prec_setup,
                                    jac_full_prec_solveV3);
    if(check_cvode_flag(&flag, "CVSpilsSetPreconditioner", 1)) exit(-1);

    flag = CVSpilsSetEpsLin(systemParam.cvodeMemPtr, epslin);
    if (check_cvode_flag(&flag, "CVSpilsSetEpsLin", 1)) exit(-1);
  }

  /* Set the maximum number of internal steps per CVode call and the maximum
   * allowable internal steps. */
  flag = CVodeSetMaxNumSteps(systemParam.cvodeMemPtr, maxsteps);
  if (check_cvode_flag(&flag, "CVodeSetMaxNumSteps", 1)) exit(-1);
  flag = CVodeSetMaxStep(systemParam.cvodeMemPtr, maxdt_internal);
  if (check_cvode_flag(&flag, "CVodeSetMaxStep", 1)) exit(-1);

  flag = CVodeSetNonlinConvCoef(systemParam.cvodeMemPtr, nlconvcoef);
  if (check_cvode_flag(&flag, "CVodeSetNonlinConvCoef", 1)) exit(-1);

  N_VDestroy_Serial(y0);

#ifdef USE_MPI
  systemParam.rank_weights =(double *)malloc(sizeof(double)*systemParam.nranks);
#ifdef ZERORK_GPU_ONLY
  double my_weight = systemParam.gpu_id >= 0 ? 1.0 : 0.0; //GPU's do all the work
  int n_gpus_global = 0;
  MPI_Allreduce(&n_gpus,&n_gpus_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if(n_gpus_global == 0) my_weight = 1.0;
#else
  double my_weight = systemParam.gpu_id >= 0 ? 10.0 : 1.0;
#endif
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
void zerork_reactor_setup_defaults
(
    double rtol,
    double atol,
    double sparsethresh
)
{
    int abstol_dens = 0; //mass fraction based abstol
    int CVCP = 1;  //Const pressure derivative (using Converge equation).
    int maxsteps = 1e6;
    int maxord = 5;  //CVODE Default
    double maxdt_internal = DBL_MAX; //TODO: Reasonable value?
    double nlconvcoef = 0.1; //CVODE Default
    double epslin = 0.05;    //CVODE Default
    int multireac = 0;  //No GPU
    int gpu_id = -1;    //No GPU
    int n_reactors_max = 1024; //Won't be used
    int n_reactors_min = 256;  //Won't be used
    int lsolver = 4; //Sparse Iterative Preconditioned

    //Won't be used (until GPU is on by default):
    double logiA = 0.0;
    double logiB = 0.0;
    double logiK = 0.0;
    double logiQ = 0.0;
    double logiM = 0.0;
    double logirNu = 0.0;

    const char * ck_log_filename = "";       //No output
    const char * reactor_log_filename = "";  //No output

    int verbosity = 0;

    int cusolver_factor_alg = 2;
    int cusolver_solve_alg = 0;

    zerork_reactor_setup_full(verbosity,maxsteps,maxord,rtol,atol,abstol_dens,
                           maxdt_internal,sparsethresh,nlconvcoef,epslin,
                           lsolver,CVCP,"mech.dat","therm.dat",ck_log_filename,
                           reactor_log_filename,n_reactors_max,n_reactors_min,
                           logiA,logiK,logiQ,logiM,logiB,logirNu,
                           cusolver_factor_alg, cusolver_solve_alg,
                           &multireac,&gpu_id);
}


#ifndef ZERORK_CONVERGE_RELEASE
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
void zerork_solve_reactors_mz(const int nReactors, const double dt,
                           const double *T, const double *P,
                           const double *dpdt, const double *volume,
                           double *massFracPtr, double *reactorCost,
                           double *reactorGpu)
{
  static int n_calls = 0;
  n_calls++;

  int nzones = 0;

  //Worst case is that every reactor is it's own zone.
  // TODO: This could be much higher in parallel, need to let zoner allocate
  // N.B. memory requirement to ensure these are large enough
  //      mean we have to allocate after counting up zones.
  std::vector<double> temp_zones;
  std::vector<double> press_zones;
  std::vector<double> dpdt_zones;
  std::vector<double> massfracs_zones;
  std::vector<double> massfracs_zones0;
  std::vector<double> cost_zones;
  std::vector<double> gpu_flag_zones;

  //Zone me up
  multizone_make_zones(*systemParam.mech, nReactors, T, P, dpdt, volume,
                       massFracPtr, reactorCost, reactorGpu, nzones,
                       temp_zones, press_zones,
                       dpdt_zones, massfracs_zones, cost_zones, gpu_flag_zones);

  massfracs_zones0 = massfracs_zones;

#ifdef USE_MPI
  //Every processor has all zones
  //  Pick zones to do and do them
  int zones_per_rank = nzones / systemParam.nranks;
  int remaining_zones = nzones % systemParam.nranks;
  std::vector<int> zone_start(systemParam.nranks+1);
  zone_start[0] = 0;
  for(int irank = 0; irank < systemParam.nranks; ++irank)
  {
    int zones_on_rank = zones_per_rank;
    if(irank < remaining_zones) zones_on_rank += 1;
    zone_start[irank+1] = zone_start[irank] + zones_on_rank;
  }
  assert(zone_start[systemParam.nranks] == nzones);
  int local_zone_count = zone_start[systemParam.rank+1] - zone_start[systemParam.rank];
#else
  int local_zone_count = nzones;
  std::vector<int> zone_start(2);
  zone_start[0] = 0;
  zone_start[1] = nzones;
#endif


  zerork_solve_reactors(local_zone_count, dt,
                     &temp_zones[zone_start[systemParam.rank]],
                     &press_zones[zone_start[systemParam.rank]],
                     &dpdt_zones[zone_start[systemParam.rank]],
                     &massfracs_zones[zone_start[systemParam.rank]*systemParam.nSpc],
                     &cost_zones[zone_start[systemParam.rank]],
                     &gpu_flag_zones[zone_start[systemParam.rank]]);

#ifdef USE_MPI
  //Broadcast zone results
  for(int irank = 0; irank < systemParam.nranks; ++irank)
  {
    int start = zone_start[irank];
    int count = zone_start[irank+1] - start;
    if(count > 0) {
      MPI_Bcast(&massfracs_zones[start*systemParam.nSpc],count*systemParam.nSpc,
                MPI_DOUBLE,irank,MPI_COMM_WORLD);
      MPI_Bcast(&cost_zones[start],count,
                MPI_DOUBLE,irank,MPI_COMM_WORLD);
    }
  }
#endif


  //Remap to reactors
  multizone_remap(*systemParam.mech,
                   nzones,
                   massfracs_zones0,
                   massfracs_zones,
                   nReactors,
                   massFracPtr,
                   cost_zones, reactorCost,
                   gpu_flag_zones, reactorGpu);
}

extern "C"
void zerork_solve_reactors(const int nReactors, const double dt,
                        const double *T, const double *P,
                        const double *dpdt,
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
  int have_gpu = systemParam.gpu_id > -1 ? 1 : 0;
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
      if(systemParam.rank==0 && systemParam.verbosity > 0) printf("Rank[%d] reactors_wanted, rank_weight: %d, %f\n",i,reactors_wanted[i],systemParam.rank_weights[i]);
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
    for(int i = 0; i < systemParam.nranks; ++i)
    { //note where we are in the list
      if(reactor_deficit[i].first == systemParam.rank)
      {
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

#ifdef ZERORK_DUMP_REACTORS

  static int nreactors_all_time = 0;
  int reactor_label = nreactors_all_time;
  for(int irank = 0; irank < systemParam.nranks; ++irank)
  {
    if(irank < systemParam.rank)
    {
      reactor_label += reactors_on_rank[irank];
    }
  }
  for(int ireactor=0; ireactor < reactors_on_rank[systemParam.rank]; ++ireactor)
  {
      reactor_label += 1;
      if(reactor_label < 1e8)
      {
          char buff[9];
          snprintf(buff,9,"%08d",reactor_label);
          std::string stateFile = std::string("zerork_state")+buff;
          FILE* dumpFile=fopen(stateFile.c_str(),"w");
          double currReactorTemp = T[ireactor];
          double currReactorPress = P[ireactor];
          for(int j = 0; j<systemParam.nSpc; ++j)
          {
            fprintf(dumpFile,"%40.32g\n",massFracPtr[ireactor*systemParam.nSpc+j]);
          }
          fprintf(dumpFile,"%40.32g\n",T[ireactor]);
          fprintf(dumpFile,"%40.32g\n",P[ireactor]);
          fprintf(dumpFile,"%40.32g\n",reactorCost[ireactor]);
          fprintf(dumpFile,"%40.32g\n",reactorGpu[ireactor]);
          fclose(dumpFile);
      }
  }

#ifdef USE_MPI
  nreactors_all_time += total_reactors;
#else
  nreactors_all_time += nreactors_calc;
#endif
#endif

  //Allocate space for reactors
  int nreactors_alloc = max(nReactors,nreactors_calc);
  std::vector<double> T_loc(nreactors_alloc);
  std::vector<double> P_loc(nreactors_alloc);
  std::vector<double> dpdt_loc(nreactors_alloc);
  std::vector<double> mf_loc(nreactors_alloc*systemParam.nSpc);
  std::vector<double> nstep_loc(nreactors_alloc);
  std::vector<double> gpu_flag_loc(nreactors_alloc);
  std::vector<double> cost_loc(nreactors_alloc);

  std::vector<double> reactorTime(nreactors_calc,0.0);
  std::vector<int> gpuSolve(nreactors_calc,0);

  //TODO: Make sorting optional and/or test effect of sorting
  bool sort_reactors = false;
  std::vector<indexed_double> dtEst;
  if(sort_reactors)
  {
    startTime = getHighResolutionTime();
    dtEst.resize(nreactors_alloc);
    for(int k = 0; k < nReactors; ++k)
    {
//    double dt_k = cvode_dti(&systemParam, T[k], P[k], dpdt[k],
//              dt, &(massFracPtr[k*systemParam.nSpc]));
//      double dt_k = chemeq2_dt(&systemParam, T[k], P[k],
//                               &(massFracPtr[k*systemParam.nSpc]));
//      dtEst[k] = std::make_pair(k,dt_k/dt); //1/(number of steps to complete)
      //TODO: Probably not the right place to put this.
      if(reactorCost[k] == 0.0) reactorCost[k] = 1.0;
      dtEst[k] = std::make_pair(k,reactorCost[k]); //1/(number of steps to complete)
    }
    //sort(&dtEst[0],&dtEst[nReactors],lesser_second<int,double>()); //Smallest first
    sort(dtEst.begin(),dtEst.end(),greater_second<int,double>()); //Biggest first
    for(int k = 0; k < nReactors; ++k)
    {
      T_loc[dtEst[k].first] = T[k];
      P_loc[dtEst[k].first] = P[k];
      dpdt_loc[dtEst[k].first] = dpdt[k];
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
    memcpy(&dpdt_loc[0],dpdt,sizeof(double)*nReactors);
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
      int recv_count_per_reactor = 6+systemParam.nSpc;
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
        memcpy(&dpdt_loc[recv_idx], &recv_buf[buf_idx],
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
      int send_count_per_reactor = 6+systemParam.nSpc;
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
            memcpy(&send_buf[buf_idx], &dpdt_loc[send_idx],
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
    std::vector<double> dpdt_loc_tmp(nreactors_calc);
    std::vector<double> nstep_loc_tmp(nreactors_calc);
    std::vector<double> gpu_flag_loc_tmp(nreactors_calc);
    std::vector<double> cost_loc_tmp(nreactors_calc);
    std::vector<double> mf_loc_tmp(nreactors_calc*systemParam.nSpc);
    for(int k = 0; k < nReactors; ++k)
    {
      T_loc_tmp[dtEst[k].first] = T[k];
      P_loc_tmp[dtEst[k].first] = P[k];
      dpdt_loc_tmp[dtEst[k].first] = dpdt[k];
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
      dpdt_loc_tmp[dtEst[k].first] = dpdt_loc[k];
      nstep_loc_tmp[dtEst[k].first] = nstep_loc[k];
      gpu_flag_loc_tmp[dtEst[k].first] = gpu_flag_loc[k];
      cost_loc_tmp[dtEst[k].first] = cost_loc[k];
      memcpy(&mf_loc_tmp[dtEst[k].first*systemParam.nSpc],
             &mf_loc[k*systemParam.nSpc],
             sizeof(double)*systemParam.nSpc);
    }
    memcpy(&T_loc[0],&T_loc_tmp[0],sizeof(double)*nreactors_calc);
    memcpy(&P_loc[0],&P_loc_tmp[0],sizeof(double)*nreactors_calc);
    memcpy(&dpdt_loc[0],&dpdt_loc_tmp[0],sizeof(double)*nreactors_calc);
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
    //No gpu assigned to this thread.  Solve in serial on cpu.
    for(int k = 0; k < nreactors_calc; ++k)
    {
        startTime=getHighResolutionTime();
        nstep_loc[k] = zerork_reactor_solve(dt, T_loc[k], P_loc[k], dpdt_loc[k],
                             &(mf_loc[k*systemParam.nSpc]));
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
#ifdef ZERORK_GPU_ONLY
      const int update_method = -1;
#else
      const int update_method = 1;
#endif
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
  //memcpy(reactorCost,&nstep_loc[0],sizeof(double)*nReactors);
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


#endif //ZERORK_CONVERGE_RELEASE


long int zerork_reactor_solve(const double dt, const double T,
                           const double P, const double dpdt,
                           double *massFracPtr)
{
  int flag;
  double tnext,tcurr;

  // reset the mass fraction, density and relative vol at new temperature
  systemParam.currReactor = 0;
  systemParam.Press[0] = P;
  systemParam.invDens[0] = systemParam.mech->getDensityFromTPY(T,P,massFracPtr);
  systemParam.invDens[0]=1.0/systemParam.invDens[0];
  systemParam.dpdt[0]=dpdt;

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
      double reactor_density = 1.0/systemParam.invDens[0];
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

  //Go for the gusto.
  //flag = CVodeSetInitStep(systemParam.cvodeMemPtr, dt);
  flag = CVodeSetMaxStep(systemParam.cvodeMemPtr,
                         std::min(systemParam.maxdt_internal,tnext));

  int maxRetry = 6;
  int nRetry   = 0;
#if 1
  while(1)
  {
      flag = CVode(systemParam.cvodeMemPtr, tnext, systemState, &tcurr, CV_NORMAL);

      if (check_cvode_flag(&flag, "CVode", 1)) {
        printf("WARNING: Failed integration current temperature: %g\n",NV_Ith_S(systemState,systemParam.nSpc));
//        for(int j = 0; j<= systemParam.nSpc; ++j)
//        {
//          printf("%g\n", NV_Ith_S(systemState,j));
//        }
//        double press = systemParam.mech->getPressureFromTVY(
//            NV_Ith_S(systemState,systemParam.nSpc),
//            systemParam.invDens[0],
//            NV_DATA_S(systemState));
//        printf("%g\n",press);
//        exit(-1);
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
#else
  while(tcurr<tnext)
  {
      flag = CVode(systemParam.cvodeMemPtr, tnext, systemState, &tcurr, CV_ONE_STEP);
      printf("%g\t",tcurr);
      for(int j = 0; j<= systemParam.nSpc; ++j)
      {
        printf("%g\t", NV_Ith_S(systemState,j));
      }
      printf("\n");

      if (check_cvode_flag(&flag, "CVode", 1)) {
        printf("WARNING: Failed integration current Temp and  %g\n",NV_Ith_S(systemState,systemParam.nSpc));
        for(int j = 0; j<= systemParam.nSpc; ++j)
        {
          printf("%g\n", NV_Ith_S(systemState,j));
        }
        double press = systemParam.mech->getPressureFromTVY(
            NV_Ith_S(systemState,systemParam.nSpc),
            systemParam.invDens[0],
            NV_DATA_S(systemState));
        printf("%g\n",press);
        exit(-1);
      }
  }
#endif
  zerork_update_cvode_stats(&systemParam, systemParam.cvodeMemPtr, false);
  systemParam.simTime += getHighResolutionTime() - startTime;
  long int nsteps;
  CVodeGetNumSteps(systemParam.cvodeMemPtr,&nsteps);

  //Failed to complete
  if(flag < 0)
  {
    printf("WARNING: Failed to complete CPU integration.\n");
    printf("         Check results carefully.\n");
  }

#ifdef ZERORK_FULL_DEBUG
//  if(flag < 0)
//  {
//    printf("Failed to complete integration. Printing state and quitting.\n");
//    printf("WARNING: Failed to complete CPU integration.\n");
//    printf("         Check results carefully.\n");
//    for(int j = 0; j<systemParam.nSpc;++j)
//    {
//      printf("%30.24g %30.24g %30.24g\n",
//             massFracPtr[j],
//             NV_Ith_S(systemState,j),
//             NV_Ith_S(systemState,j) - massFracPtr[j]);
//    }
//    printf("T: %30.24g\n",T);
//    printf("P: %30.24g\n",P);
//    exit(-1);
//  }
#endif

  memcpy(massFracPtr,NV_DATA_S(systemState),sizeof(double)*systemParam.nSpc);

  N_VDestroy_Serial(systemState);
  return nsteps;
}


extern "C" void zerork_species_reaction_rates(const double T, const double P,
                                const double * molar_concentration_cgs,
                                double * wdot)
{
  double * molar_concentration_mks = (double*)malloc(sizeof(double)*systemParam.nSpc);
  for(int j=0; j<systemParam.nSpc; j++)
  {
    molar_concentration_mks[j] = 1000.0*molar_concentration_cgs[j];
  }

  // compute the molar production rates at the current state (aka wdot)
  systemParam.mech->getReactionRates(T,molar_concentration_mks,
                                     wdot,systemParam.createRate,
                                     systemParam.destroyRate,systemParam.fwdROP);

  // wdot = [kmol/m^3/s] * (1000 mol/ 1 kmol) * (1 m / cm 100)^3 = [mol/cm^3/s]
  for(int j=0; j<systemParam.nSpc; j++)
  {
    wdot[j] /= 1000.0;
  }
  free(molar_concentration_mks);
}

#ifdef ZERORK_USE_SAGE_RHS
extern "C"
int sage_fun(realtype t, N_Vector y, N_Vector ydot, void *f_data);
#include "../../zerork_converge_udf/const_shared.h" //for PRECISION and real
#include "../../zerork_converge_udf/structures_shared.h" //for mechanidsm_definition
extern "C"
struct mechanism_definition* get_sage_mech_ptr();

int const_dpdt_wsr_sage(realtype t, N_Vector y, N_Vector ydot,
			 void *user_data)
{
  cv_param *cvp=(cv_param *)user_data;
  double * y_ptr = NV_DATA_S(y);
  double * ydot_ptr = NV_DATA_S(ydot);
  double startTime=getHighResolutionTime();

  double Temp = y_ptr[cvp->nSpc]*cvp->Tref;
  double currPress, invDensCurr;
  //ASSUME : t0 == 0.0;
  if(systemParam.constPress)
  {
#ifdef CONVERGE_CONST_PRESS
    currPress = cvp->Press[cvp->currReactor];
#else
#ifdef CONVERGE_DPDT_UNITS
    currPress = cvp->Press[cvp->currReactor]+(t)*cvp->dpdt[cvp->currReactor]/10.0;
#else
    currPress = cvp->Press[cvp->currReactor]+(t)*cvp->dpdt[cvp->currReactor];
#endif
#endif
    invDensCurr = 1.0/(cvp->mech->getDensityFromTPY(Temp,currPress,y_ptr));
  }
  else
  {
    invDensCurr = cvp->invDens[cvp->currReactor];
    currPress = cvp->Press[cvp->currReactor]; //TODO: FIXME
  }


  // set concentration via density and mass fraction
  cvp->mech->getCfromVY(invDensCurr,y_ptr,cvp->conc);
  N_Vector sage_state = N_VNew_Serial(cvp->nSpc+1); // molar concentration in cgs and Temp
  N_Vector sage_deriv = N_VNew_Serial(cvp->nSpc+1); //molar production rate in cgs
  double * sage_state_ptr = NV_DATA_S(sage_state);
  double * sage_deriv_ptr = NV_DATA_S(sage_deriv);
  for(int j=0; j<cvp->nSpc; ++j)
  {
    sage_state_ptr[j] = cvp->conc[j]/1000.0; //mks -> cgs
  }
  sage_state_ptr[cvp->nSpc] = Temp;
  struct mechanism_definition* sage_mech=get_sage_mech_ptr();

  //Set the state for sage
  sage_mech->temp = sage_state_ptr[cvp->nSpc];
  sage_mech->pres = currPress*10.0; //mks -> cgs
  sage_mech->dpdt = cvp->dpdt[cvp->currReactor]*10.0; //mks -> cgs
#ifdef CONVERGE_DPDT_UNITS
  sage_mech->dpdt *= 10.0;
#endif
  sage_mech->solve_temp = 1;

  if(systemParam.constPress)
  {
    sage_mech->conv_flag = 0;
    sage_mech->conp_flag = 1;
  } else {
    sage_mech->conv_flag = 1;
    sage_mech->conp_flag = 0;
  }

  sage_fun(t,sage_state,sage_deriv,sage_mech);


  for(int j=0; j<cvp->nSpc; ++j)
  {
    // cgs molar production -> mks -> mass fraction derivative
    ydot_ptr[j]=(sage_deriv_ptr[j]*1000.0)*(cvp->molWt[j])*(invDensCurr);
  }
  ydot_ptr[cvp->nSpc]=sage_deriv_ptr[cvp->nSpc]/cvp->Tref;


  N_VDestroy_Serial(sage_state);
  N_VDestroy_Serial(sage_deriv);

  (cvp->nFunc)++;
  cvp->funcTime += getHighResolutionTime() - startTime;
  return 0;
}
#endif //ZERORK_USE_SAGE_RHS


extern "C" void zerork_reactor_close(void)
{
  // Free integrator memory
  CVodeFree(&(systemParam.cvodeMemPtr));

  // Free parameter memory
  free(systemParam.Press);
  free(systemParam.invDens);
  free(systemParam.meanCvMass);
  free(systemParam.dpdt);
  free(systemParam.dTemp_dt);
  free(systemParam.molWt);
  free(systemParam.invMolWt);
  free(systemParam.netProd);
  free(systemParam.fwdROP);
  free(systemParam.createRate);
  free(systemParam.destroyRate);
  free(systemParam.conc);
  free(systemParam.Energy);
  free(systemParam.CvMass);


  {
    free(systemParam.systemState_host);
    free(systemParam.systemDeriv_host);
    free(systemParam.tmp1_host);
    free(systemParam.tmp2_host);
    free(systemParam.tmp3_host);
  }



  free_Jsparse(systemParam.sparseMtx);
  delete systemParam.mech;

#ifdef USE_MPI
  free(systemParam.rank_weights);
#endif

  // clean up procedures
  fclose(reactor_log_file);

//  printf("Exiting zerork_reactor_close()\n");
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

static double calc_max_dt_ratio(double dt)
{
    //double max_dt_diff_log = 2.20/(1+3.2*exp(-log_dt_first-2.2))+0.015;
    return
        systemParam.logiA + //Minimum of logistic function
        (systemParam.logiK -systemParam.logiA) /   //Maximum of logistic func.
        (
            pow((1.0+systemParam.logiQ* //Pre-exponential
                 exp(-systemParam.logiB*(dt-systemParam.logiM)) //Rate of change and zero crossing
            ),systemParam.logirNu) //Growth rate bias
        );
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
	     systemParam.sparseMtx->reduceNNZ[0],
	     systemParam.sparseMtx->LUnnz[0]);
  printf("  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %8.4e\n",
             systemParam.jacSetupTime,
             systemParam.precSetupTime,
             systemParam.colPermTime,
             systemParam.jacFactorTime,
             systemParam.backsolveTime,
             systemParam.funcTime, otherTime, systemParam.simTime);
  ++print_count;
}

/*
static void zerork_update_cvode_stats(cv_param* cvp, void* cvode_mem_ptr)
{
  int flag;
  long int nsteps,nfevals,nlinsetups,netfails;
  long int nniters,nncfails,ngevals;
  int qlast, qcur;
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
*/

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

  if(systemParam.lsolver < 2)
  {
    double cv_setup_time = 0.0;
    double cv_solve_time = 0.0;
    double cv_jac_eval_time = 0.0;
    {
      //Note we reset the mr pointer every time but not the single reactor
      CVUserTimedGetSetupTime(cvode_mem_ptr, &cv_setup_time);
      CVUserTimedGetSolveTime(cvode_mem_ptr, &cv_solve_time);
      CVUserTimedGetJacEvalTime(cvode_mem_ptr, &cv_jac_eval_time);
      CVDlsGetNumRhsEvals(cvode_mem_ptr, &nfeDQ);
      cvp->jacSetupTime  = cv_jac_eval_time; //N.B. disagree-able nomenclature
      cvp->jacFactorTime = cv_setup_time;
      cvp->backsolveTime = cv_solve_time;
    }
  }
  //TODO: HACK
  if(systemParam.lsolver == 1)
  {
    cvp->funcTime *=  (cvp->nfevals)/((double)(cvp->nfevals+nfeDQ));
  }
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

