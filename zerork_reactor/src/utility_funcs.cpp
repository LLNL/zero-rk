#include "utility_funcs.h"

#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include "cv_param_sparse.h"


double get_val_sp_matrix(int* colSum, int* rowIdx, double* vals, int i, int j)
{
  int n;
  for(n = colSum[j]; n < colSum[j+1]; ++n)
  {
      if(rowIdx[n] == i) {return vals[n];}
      else if(rowIdx[n] > i) {return 0.0;}
  }
  return 0.0;
}

double get_val_sp_matrix_csr(int* rowSum, int* colIdx, double* vals, int i, int j)
{
  int n;
  for(n = rowSum[i]; n < rowSum[i+1]; ++n)
  {
      if(colIdx[n] == j) {return vals[n];}
      else if(colIdx[n] > j) {return 0.0;}
  }
  return 0.0;
}


void print_sp_matrix(int m, int n, int* colSum, int* rowIdx, double* avals)
{
    int i,j;
    for(i = 0; i<m; ++i)
    {
        for(j = 0; j < n; ++j)
        {
           double val = get_val_sp_matrix(colSum,rowIdx,avals,i,j);
           printf("%14.8e\t",val); 
        }
        printf("\n");
    }
}


void print_sp_matrix_csr(int m, int n, int* rowSum, int* colIdx, double* avals)
{
    int i,j;
    for(i = 0; i<m; ++i)
    {
        for(j = 0; j < n; ++j)
        {
           double val = get_val_sp_matrix_csr(rowSum,colIdx,avals,i,j);
//           printf("%14.8e\t",val); 
           printf("%10.3e\t",val); 
        }
        printf("\n");
    }
}



void permute_sparse_csc(int n,const int* aColSum, const int *aRowIdx, const double *aVals, int* bColSum, int* bRowIdx, double *bVals, int* perm_c)
{
// * perm_c: Column permutation vector of size n, which defines the 
// *         permutation matrix Pc; perm_c[i] = j means column i of A is 
// *         in position j in A*Pc.
  int j,k,acol,bnnz;
  bnnz = 0;
  bColSum[0] = bnnz;
  for(j = 0; j < n; ++j) //columns of permuted matrix
    {
      bColSum[j+1] = bColSum[j]; //start the counter
      acol = perm_c[j]; //column of original matrix
      for(k = aColSum[acol]; k < aColSum[acol+1]; ++k)
        {
          bRowIdx[bnnz] = aRowIdx[k];
          bVals[bnnz] = aVals[k];
          ++bnnz;
          bColSum[j+1] = bnnz;
        }
    }
}


void permute_sparse_csc_iperm(int n,const int* aColSum, const int *aRowIdx, const double *aVals, int* bColSum, int* bRowIdx, double *bVals, int* iperm_c)
{
  int j;
  int* perm_c = new int[n];
  for(j = 0; j < n; ++j)
  {
    perm_c[iperm_c[j]] = j;
  }
  permute_sparse_csc(n, aColSum, aRowIdx, aVals, bColSum, bRowIdx, bVals, perm_c);
  delete [] perm_c;
}


double getHighResolutionTime(void)
{
    struct timeval tod;

    gettimeofday(&tod, NULL);
    double time_seconds = (double) tod.tv_sec + ((double) tod.tv_usec / 1000000.0);
    return time_seconds;
}


/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */
int check_cvode_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}



double chemeq2_dt(void *user_data, double T,
                  double pres, double * massFrac)
{
  zerork::mechanism* mech = ((cv_param *) user_data)->mech;
  int nSpc = mech->getNumSpecies();
  int nStep = mech->getNumSteps();

  double invDens = 1.0/mech->getDensityFromTPY(T, pres, massFrac);
  double * conc = new double[nSpc];
  double * netProd = new double[nSpc];
  double * createRate = new double[nSpc];
  double * destroyRate = new double[nSpc];
  double * fwdROP = new double[nStep];
  double * molWt = new double[nSpc];
  double * ydot = new double[nSpc];
  double * q = new double[nSpc];
  double * p = new double[nSpc];
  mech->getCfromVY(invDens,massFrac,conc);
  mech->getReactionRates(T,conc,netProd,createRate,
			      destroyRate,fwdROP);
  mech->getMolWtSpc(molWt);

  double epsilon = 1.0e-8;
  double alpha = 0.5*sqrt(epsilon);
  double dtmin = 1.0e+300;
  double dtj   = 1.0e+300;
  for(int j = 0; j < nSpc; ++j)
  {
    ydot[j]=netProd[j]*molWt[j]*invDens;
    q[j]=createRate[j]*molWt[j]*invDens;
    p[j]=(q[j] - ydot[j])/massFrac[j];
    double c1 = 0.1*epsilon*q[j];
    if( c1 > p[j]*massFrac[j] )
    {
      dtj = alpha/(p[j] + 1.0e-300);
    }
    else
    {
      dtj = (alpha*massFrac[j])/(fabs(ydot[j])+1.0e-300);
      if(dtj == 0.0) dtj = 1.0e+300;
    }
    dtmin = min(dtj,dtmin);
  }

  delete [] conc;
  delete [] netProd;
  delete [] createRate;
  delete [] destroyRate;
  delete [] fwdROP;
  delete [] molWt;
  delete [] ydot;
  delete [] q;
  delete [] p;

  return dtmin;
}

#include "cvode/cvode.h"
#include "cvode/cvode_dense.h"
#include "cvode/cvode_spgmr.h"
#include <nvector/nvector_serial.h>
#include "ode_funcs.h"
#include "matrix_funcs.h"

double cvode_dti(void *user_data, double T, double pres,
                 double dpdt, double dt, double *massFrac)
{
  cv_param *cvp = (cv_param *)user_data;
  int flag;
  double tnext,tcurr;

  // reset the mass fraction, density and relative vol at new temperature
  int currReactor_old = cvp->currReactor; //save old status of currReactor

  cvp->currReactor = 0; //do all ops on reactor 0
  double P_old,invDens_old,dpdt_old;
  P_old = cvp->Press[cvp->currReactor];
  invDens_old = cvp->invDens[cvp->currReactor];
  dpdt_old = cvp->dpdt[cvp->currReactor];

  cvp->Press[cvp->currReactor] = pres;
  cvp->invDens[cvp->currReactor] = 1.0/cvp->mech->getDensityFromTPY(T,pres,massFrac);
  cvp->dpdt[cvp->currReactor]=dpdt;

  //TODO: can we eliminate this memcpy
  N_Vector systemState = N_VNew_Serial(cvp->nSpc+1);
  memcpy(NV_DATA_S(systemState),massFrac,sizeof(double)*cvp->nSpc);
  NV_Ith_S(systemState,cvp->nSpc) = T/cvp->Tref;

  // reset the time
  tcurr=0.0;
  tnext=dt;

  void* cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  flag = CVodeInit(cvode_mem, const_vol_wsr, tcurr, systemState);
  flag = CVodeSetErrFile(cvode_mem, NULL); //Turn off error printing
                                           //since we expect to fail
                                           //in matrix setup
  flag = CVodeSStolerances(cvode_mem,cvp->rtol,cvp->atol);
  flag = CVodeSetUserData(cvode_mem, cvp);

#if 0
    flag = CVDense(cvode_mem, cvp->nSpc+1);
    if(check_cvode_flag(&flag, "CVDense", 1)) exit(-1);
#else
    flag = CVSpgmr(cvode_mem, PREC_LEFT, 5);//krylov_dim
    if(check_cvode_flag(&flag, "CVSpgmr", 1)) exit(-1);

    flag = CVSpilsSetGSType(cvode_mem, MODIFIED_GS);
    if(check_cvode_flag(&flag, "CVSpilsSetGSType", 1)) exit(-1);

    flag = CVSpilsSetPreconditioner(cvode_mem, jac_full_prec_setup,
    				  jac_full_prec_solveV3);
    if(check_cvode_flag(&flag, "CVSpilsSetPreconditioner", 1)) exit(-1);
#endif

  cvp->sparseMtx->abort_matrix_eval=true;
  flag = CVode(cvode_mem,tnext,systemState,&tcurr,CV_ONE_STEP);
  cvp->sparseMtx->abort_matrix_eval=false;
  if (flag != CV_SUCCESS && flag != CV_LSETUP_FAIL)
  {
    check_cvode_flag(&flag, "CVode", 1);
  }

  //Get the step cvode would have taken
  double dt_init;
  flag = CVodeGetActualInitStep(cvode_mem,&dt_init);

  CVodeFree(&cvode_mem);
  N_VDestroy_Serial(systemState);

  cvp->Press[cvp->currReactor] = P_old;
  cvp->invDens[cvp->currReactor] = invDens_old;
  cvp->dpdt[cvp->currReactor] = dpdt_old;
  cvp->currReactor = currReactor_old;

  return dt_init;
}


