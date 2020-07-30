#include <limits>

#include "matrix_funcs.h"
#include "cv_param_sparse.h"
#include "ode_funcs.h"
#include "utility_funcs.h"

#if defined  SUNDIALS2
int jac_full_prec_setup(realtype t, N_Vector y, N_Vector fy,
			booleantype jok, booleantype *jcurPtr,
			realtype gamma, void *user_data,
			N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
#elif defined SUNDIALS3 || defined SUNDIALS4
int jac_full_prec_setup(realtype t, N_Vector y, N_Vector fy,
			booleantype jok, booleantype *jcurPtr,
			realtype gamma, void *user_data)
#endif
{
  cv_param *cvp=(cv_param *)user_data;
  double startTime;
  int lwork=0;      // SuperLU allocates its own memeory
  void *work=NULL;  // ''
  double rpg,rcond; // reciprocal pivot growth, recip. condition number
  double *ferr = NULL; //Not factoring so un-used
  double *berr = NULL;

  static int lunnz_last_perm = NV_LENGTH_S(y)*NV_LENGTH_S(y);
  static double gamma_last_perm = 1.0; //Won't be used until jok
  static double ff_last_perm = 100.0;
  static int num_calls = 0;

  ++num_calls;

  int misMatch=0;
  int j,flag;
  if(cvp->sparseMtx->abort_matrix_eval) { return -1; }

  if(cvp->sparseMtx->isFirstFactor) {
    lunnz_last_perm = NV_LENGTH_S(y)*NV_LENGTH_S(y);
    ff_last_perm = 100.0;
    gamma_last_perm = 1.0;
  }

  if(!jok) {
    (*jcurPtr)=TRUE; //indicate that Jacobian data was recomputed
#if defined SUNDIALS2
    setupJacobianSparse(t,y,fy,cvp,tmp1,tmp2,tmp3);
#elif defined SUNDIALS3 || defined SUNDIALS4
    int nsize = cvp->sparseMtx->nSize;
    N_Vector tmp1, tmp2, tmp3;
    tmp1 = N_VNew_Serial(nsize);
    tmp2 = N_VNew_Serial(nsize);
    tmp3 = N_VNew_Serial(nsize);
    setupJacobianSparse(t,y,fy,cvp,tmp1,tmp2,tmp3);
    N_VDestroy_Serial(tmp1);
    N_VDestroy_Serial(tmp2);
    N_VDestroy_Serial(tmp3);
#endif
    // At this point the entire Jacobian matrix should be constructed
  } else {
    (*jcurPtr)=FALSE; //indicate that Jacobian data was not recomputed
    if(cvp->sparseMtx->fakeUpdate) {return 0;}
  }

  startTime = getHighResolutionTime();
  bool bigOrderChange = fabs(log10(gamma/gamma_last_perm)) > cvp->sparseMtx->maxGammaChangeOrder;
  bool strictSamePattern = cvp->sparseMtx->strictSamePattern && jok && !bigOrderChange;
  bool doPatternCheck = (!(cvp->sparseMtx->isFirstFactor) && cvp->sparseMtx->permThresh > 0.0);
  misMatch = sparseUpdateAndCheckPattern(cvp->sparseMtx->offDiagThreshold,
      		    cvp->sparseMtx->nSize,
      		    cvp->sparseMtx->nNonZero,
      		    cvp->sparseMtx->mtxData,
      		    cvp->sparseMtx->mtxRowIdx,
      		    cvp->sparseMtx->mtxColSum,
      		    &cvp->sparseMtx->reduceNNZ,
      		    cvp->sparseMtx->reduceData,
      		    cvp->sparseMtx->reduceRowIdx,
      		    cvp->sparseMtx->reduceColSum,
                    gamma,
                    doPatternCheck,
                    strictSamePattern);

    cvp->sparseMtx->optionSLU.Fact=DOFACT;
    if(
        strictSamePattern 
        || (
              misMatch < (int)((double)cvp->sparseMtx->nSize * cvp->sparseMtx->permThresh)
              && cvp->sparseMtx->numPermReUses < 40
              && !(cvp->sparseMtx->isFirstFactor)
              && !(bigOrderChange)
           )
      )
    {
         cvp->sparseMtx->optionSLU.Fact=SamePattern;
         if(!strictSamePattern && jok) {++(cvp->sparseMtx->numPermReUses);}
         ++(cvp->sparseMtx->numPermReUses);
    } else {
        cvp->sparseMtx->numPermReUses = 0;
    }

    cvp->precSetupTime += getHighResolutionTime() - startTime;


    if(!(cvp->sparseMtx->isFirstFactor)  && cvp->sparseMtx->optionSLU.Fact != SamePattern_SameRowPerm)
    {
            Destroy_SuperNode_Matrix(&(cvp->sparseMtx->Lslu));
            Destroy_CompCol_Matrix(&(cvp->sparseMtx->Uslu));
    }

    //Update number of non-zeros in M.
    ((NCformat *)cvp->sparseMtx->Mslu.Store)->nnz = cvp->sparseMtx->reduceNNZ;

    cvp->sparseMtx->Bslu.ncol=0; // in dlinsolx1.c example this is supposed to
                                 // indicate that a solution is not needed only
                                 // the factorization

    // Perform permutation
    if(cvp->sparseMtx->optionSLU.Fact == DOFACT)
    {
      startTime = getHighResolutionTime();
      get_perm_c(MMD_AT_PLUS_A, &(cvp->sparseMtx->Mslu), cvp->sparseMtx->colPermutation);
      ++(cvp->nColPerm); 
      cvp->colPermTime += getHighResolutionTime() - startTime;
    }

    startTime = getHighResolutionTime();
    dgssvx(&(cvp->sparseMtx->optionSLU),&(cvp->sparseMtx->Mslu),
        cvp->sparseMtx->colPermutation,cvp->sparseMtx->rowPermutation,
        cvp->sparseMtx->colElimTree,cvp->sparseMtx->equed,cvp->sparseMtx->Rvec,
        cvp->sparseMtx->Cvec,&(cvp->sparseMtx->Lslu),&(cvp->sparseMtx->Uslu),
        work,lwork,&(cvp->sparseMtx->Bslu),&(cvp->sparseMtx->Xslu),&rpg,&rcond,
        ferr,berr,
#if SUPERLU_MAJOR_VERSION > 4
        &cvp->sparseMtx->Glu,
#endif
        &(cvp->sparseMtx->mem_usage),&(cvp->sparseMtx->statSLU), &flag);
    ++(cvp->nJacFactor);

    cvp->sparseMtx->LUnnz = ((SCformat *)cvp->sparseMtx->Lslu.Store)->nnz + ((NCformat *)cvp->sparseMtx->Uslu.Store)->nnz;
    cvp->sparseMtx->fillFactor = (double)cvp->sparseMtx->LUnnz/(double)cvp->sparseMtx->reduceNNZ;

    bool redoFact = false;

    if(cvp->sparseMtx->optionSLU.DiagPivotThresh > 0.0) {
      for( j = 0; j < cvp->sparseMtx->nSize; ++j ) {
        if( cvp->sparseMtx->rowPermutation[j] == -1 ) {
          int colElems = cvp->sparseMtx->reduceColSum[j+1] - cvp->sparseMtx->reduceColSum[j];
          redoFact = true;
        }
      }
    }

    if( redoFact ) 
    {
      cvp->sparseMtx->optionSLU.Fact = DOFACT;
      Destroy_SuperNode_Matrix(&(cvp->sparseMtx->Lslu));
      Destroy_CompCol_Matrix(&(cvp->sparseMtx->Uslu));
      get_perm_c(MMD_AT_PLUS_A, &(cvp->sparseMtx->Mslu), cvp->sparseMtx->colPermutation);
      dgssvx(&(cvp->sparseMtx->optionSLU),&(cvp->sparseMtx->Mslu),
          cvp->sparseMtx->colPermutation,cvp->sparseMtx->rowPermutation,
          cvp->sparseMtx->colElimTree,cvp->sparseMtx->equed,cvp->sparseMtx->Rvec,
          cvp->sparseMtx->Cvec,&(cvp->sparseMtx->Lslu),&(cvp->sparseMtx->Uslu),
          work,lwork,&(cvp->sparseMtx->Bslu),&(cvp->sparseMtx->Xslu),&rpg,&rcond,
          ferr,berr,
#if SUPERLU_MAJOR_VERSION > 4
          &cvp->sparseMtx->Glu,
#endif
          &(cvp->sparseMtx->mem_usage),&(cvp->sparseMtx->statSLU), &flag);
      cvp->sparseMtx->LUnnz = ((SCformat *)cvp->sparseMtx->Lslu.Store)->nnz + ((NCformat *)cvp->sparseMtx->Uslu.Store)->nnz;
      cvp->sparseMtx->fillFactor = (double)cvp->sparseMtx->LUnnz/(double)cvp->sparseMtx->reduceNNZ;
    }

    cvp->jacFactorTime += getHighResolutionTime() - startTime;

    if(cvp->sparseMtx->optionSLU.Fact == DOFACT) {
      ff_last_perm = cvp->sparseMtx->fillFactor;
      lunnz_last_perm = cvp->sparseMtx->LUnnz;
      gamma_last_perm = gamma;
    }

    if(cvp->sparseMtx->isFirstFactor==1) {
      cvp->sparseMtx->isFirstFactor=0;
    }

  return flag;
}


// Solve Pz = r where P in this case is the full numerical Jacobian
//   M=I-gamma*J
#if defined SUNDIALS2
int jac_full_prec_solveV3(realtype t, N_Vector y, N_Vector fy,
		      N_Vector r, N_Vector z, realtype gamma,
		      realtype delta, int lr, void *user_data,
		      N_Vector tmp)
#elif defined SUNDIALS3 || defined SUNDIALS4
int jac_full_prec_solveV3(realtype t, N_Vector y, N_Vector fy,
                          N_Vector r, N_Vector z, realtype gamma,
                          realtype delta, int lr, void *user_data)
#endif
{
  cv_param *cvp=(cv_param *)user_data;
  int flag;
  int nsize=(cvp->nSpc)+1;
  int cpysize=nsize*sizeof(double);
  double *zptr=NV_DATA_S(z);
  double *rptr=NV_DATA_S(r);
  double startTime;

  int lwork=0;
  void *work=NULL;
  double rpg,rcond;
  double ferr[1],berr[1]; // length is the number of RHS

  cvp->sparseMtx->optionSLU.Fact=FACTORED;
  cvp->sparseMtx->Bslu.ncol=1; // in dlinsolx1.c example this is reset to the number
                               // of right hand sides

  startTime = getHighResolutionTime();
  // copy the data to the sparse matrix RHS
  memcpy(cvp->sparseMtx->vecRHS,rptr,cpysize);

  //backsolve with expert driver function
  dgssvx(&(cvp->sparseMtx->optionSLU),&(cvp->sparseMtx->Mslu),
	 cvp->sparseMtx->colPermutation,cvp->sparseMtx->rowPermutation,
	 cvp->sparseMtx->colElimTree,cvp->sparseMtx->equed,cvp->sparseMtx->Rvec,
	 cvp->sparseMtx->Cvec,&(cvp->sparseMtx->Lslu),&(cvp->sparseMtx->Uslu),
	 work,lwork,&(cvp->sparseMtx->Bslu),&(cvp->sparseMtx->Xslu),&rpg,&rcond,
	 ferr,berr,
#if SUPERLU_MAJOR_VERSION > 4
         &cvp->sparseMtx->Glu,
#endif
         &(cvp->sparseMtx->mem_usage),&(cvp->sparseMtx->statSLU), &flag);

  //copy solution to output
  memcpy(zptr,cvp->sparseMtx->vecSoln,cpysize);

  //book keeping
  ++(cvp->nBackSolve);
  cvp->backsolveTime += getHighResolutionTime() - startTime;
  return flag;
}

void setupJacobianSparse(realtype t,N_Vector y,N_Vector fy,void* user_data,N_Vector tmp1,N_Vector tmp2,N_Vector tmp3)
{
  cv_param* cvp = (cv_param *)user_data;
  double dTemp,RuTemp,multFact,startTime,startTimeDeriv,deriv_time;
  int nSpc;
  int rowId,sparseId,sparseId2,nElem;
  int j,k,N;
  double *tmp1_ptr = NV_DATA_S(tmp1);
  double *sMptr=cvp->sparseMtx->mtxData;

  startTime = getHighResolutionTime();

  nSpc=cvp->mech->getNumSpecies();
  N=nSpc+1;
  RuTemp=(cvp->mech->getGasConstant())*NV_Ith_S(y,nSpc);

  // set tmp1 to the strictly positive mass fraction array
  // set tmp2 to 1/C where C is the concentration for the strictly positive
  // mass fraction array
  for(j=0; j<nSpc; j++)
    {
       NV_Ith_S(tmp1,j)=NV_Ith_S(y,j);
       if(NV_Ith_S(tmp1,j) < cvp->minMassFrac) {
         NV_Ith_S(tmp1,j)=cvp->minMassFrac;
       }
       NV_Ith_S(tmp2,j)=cvp->molWt[j]/(NV_Ith_S(tmp1,j))*cvp->invDens;
    }

    NV_Ith_S(tmp1,nSpc)=NV_Ith_S(y,nSpc);

    // calculate the reaction info for the strictly positive case
    startTimeDeriv = getHighResolutionTime();
    const_vol_wsr(t,tmp1,tmp3,user_data);
    deriv_time = getHighResolutionTime() - startTimeDeriv;

    // SPARSE:
    // the entire sparse matrix array will be set to zero upon entry
    calcReaction_Jsparse(cvp->sparseMtx,NV_DATA_S(tmp2),cvp->fwdROP);
    // At this point sMptr stores d(wdot[k])/dC[j] ignoring the contribution
    // of perturbations in the third body species


    // ---------------------------------------------------------------------
    // compute d(Tdot)/dy[j]

    // SPARSE:
    // step 1: compute matrix vector product (d(wdot[k])/dC[j])*E[k]
    for(j=0; j<nSpc; j++) // column number
      {
        sparseId=cvp->sparseMtx->lastRowIdx[j];
        sMptr[sparseId]=0.0;
        nElem=cvp->sparseMtx->mtxColSum[j+1]-cvp->sparseMtx->mtxColSum[j];
        nElem--; // decrement the element by one since the last row is 
                 // marked
        sparseId2=cvp->sparseMtx->mtxColSum[j];
        for(k=0; k<nElem; k++)
          {
            rowId=cvp->sparseMtx->mtxRowIdx[sparseId2];
	    //printf(" k = %d, rowId = %d, sparseId = %d, sparseId2 = %d\n",k,rowId,sparseId,sparseId2); fflush(stdout);
	    sMptr[sparseId]+=sMptr[sparseId2]*(cvp->Energy[rowId]);
	    sparseId2++;
	  }
	sMptr[sparseId]*=RuTemp;
	//sMptr[sparseId]*=RuTemp; for E_RT
      }

    // step 2: add the Tdot*Cv[j] term to d(Tdot)/dy[j]
    for(j=0; j<nSpc; j++)
      {
	sparseId=cvp->sparseMtx->lastRowIdx[j];
	sMptr[sparseId]+=(cvp->dTemp_dt)*(cvp->CvMass[j]*cvp->molWt[j]);
      }

    // step 3: divide by -1/(molWt[j]*meanCvMass)
    multFact=-1.0/(cvp->meanCvMass);
    for(j=0; j<nSpc; j++)
      {
        sparseId=cvp->sparseMtx->lastRowIdx[j];
        sMptr[sparseId]*=(cvp->invMolWt[j])*multFact;
      }

    // At this point Mmtx stores d(Tdot[k])/dy[j] ignoring the contribution
    // of perturbations in the third body species

    sparseId=0;
    for(j=0; j<nSpc; j++) // jth column
      {
        nElem=cvp->sparseMtx->mtxColSum[j+1]-
              cvp->sparseMtx->mtxColSum[j];
        for(k=0; k<nElem; k++) 
          {
            rowId=cvp->sparseMtx->mtxRowIdx[sparseId];
            if(rowId<nSpc) // rowId==nSpc is the last row
              {
	        cvp->sparseMtx->mtxData[sparseId]*=
	          (cvp->molWt[rowId])*(cvp->invMolWt[j]);
	      }
	    sparseId++;
	  }
       }

    // --------------------------------------------------------------------- 
    // calculate d(ydot[k])/dT

    // step 1: perturb the temperature
    dTemp=NV_Ith_S(y,nSpc)*(cvp->sqrtUnitRnd);
    NV_Ith_S(tmp1,nSpc)+=dTemp;
    dTemp=NV_Ith_S(tmp1,nSpc)-NV_Ith_S(y,nSpc);
    multFact=1.0/dTemp;

    for(j=0; j<nSpc; j++)
      {NV_Ith_S(tmp1,j)=NV_Ith_S(y,j);}

    // step 2: calculate ydot at Temp+dTemp
    startTimeDeriv = getHighResolutionTime();
    const_vol_wsr(t,tmp1,tmp3,user_data);
    deriv_time += getHighResolutionTime() - startTimeDeriv;

    // step 3: approximate d(ydot[k])/dT with finite difference
    sparseId=cvp->sparseMtx->mtxColSum[nSpc];
    for(k=0; k<N; k++)
      {
        sMptr[sparseId]=(NV_Ith_S(tmp3,k)-NV_Ith_S(fy,k))*multFact;
        sparseId++;
      }

    cvp->jacSetupTime += getHighResolutionTime() - startTime - deriv_time;
    ++(cvp->nJacSetup);
}


void sparseOffDiagThreshCopyGamma(const double tol, const int nsize, const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int *nnzB, double B[], 
			     int BrowId[], int BcolSum[], const double gamma)
{
  int j,k,nelem;
  int AsparseId=0;
  int BsparseId=0;
  
  AsparseId=BsparseId=0;
  BcolSum[0]=0;
  for(j=0; j<nsize; j++) //
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      BcolSum[j+1]=BcolSum[j]; 
      for(k=0; k<nelem; k++)
	{
          bool diag = ArowId[AsparseId] == j;
          double currVal = -gamma*A[AsparseId];
          if (diag) currVal += 1.0;
	  if(diag || fabs(currVal) > tol) 
	    { // directly copy diagonal and any off-diagonal terms larger than
	      // tol
	      B[BsparseId]=currVal;
	      BrowId[BsparseId]=ArowId[AsparseId];
	      BcolSum[j+1]++; // increment the column sum
	      BsparseId++;    // increment B sparse index 
	    }
	  AsparseId++; // increment the A sparse index
	}
    }
  (*nnzB)=BcolSum[nsize];
}





int compare_double(const void *a, const void *b)
{
  double aval = *((double *)a);
  double bval = *((double *)b);
  return ((aval > bval) ? -1 : 1);
}

void sparseUpdateGamma(const int nsize, const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int nnzB, double B[], 
			     int BrowId[], int BcolSum[], const double gamma)
{
  int k,brow,bcol,arow,BsparseId,AsparseId,nElemAcol,lastAmatch;

  bcol = 0; 
  lastAmatch = 0;
  nElemAcol = AcolSum[bcol+1] - AcolSum[bcol];
  for(BsparseId=0; BsparseId < nnzB; ++BsparseId)
    {
       brow = BrowId[BsparseId];
       if(BsparseId >= BcolSum[bcol+1])
         {
           ++bcol;
           lastAmatch = 0;
           nElemAcol = AcolSum[bcol+1] - AcolSum[bcol];
         }
       for(k=lastAmatch; k<nElemAcol; ++k)
        {
          AsparseId = AcolSum[bcol] + k;
          arow = ArowId[AsparseId];
          if(brow == arow)
            {
              bool diag = brow == bcol;
	      B[BsparseId] = -gamma*A[AsparseId];
              if (diag) B[BsparseId] += 1.0;
              lastAmatch = k;
              break;
            }
         }
     }
}

int checkPattern(const int nsize,
                 const int nnzB1, //const double B1[],
                 const int B1rowId[], const int B1colSum[],
                 const int nnzB2, //const double B2[], 
		 const int B2rowId[], const int B2colSum[])
{
  int j,k,mismatch;
  int B1sparseId=0;
  int B2sparseId=0;

  mismatch = 0;
  for(j = 0; j<nsize; ++j) //iterating over columns
  {
     B1sparseId = B1colSum[j];
     B2sparseId = B2colSum[j];
     for(k = 0; k < nsize; ++k) //iterating over rows
     {
        if(B1sparseId == B1colSum[j+1] || B2sparseId == B2colSum[j+1]) //end of either column
        {
            mismatch += (B1colSum[j+1] - B1sparseId) + (B2colSum[j+1] - B2sparseId); //remaining terms in other matrix column
            break;
        }
        if(B1rowId[B1sparseId] < B2rowId[B2sparseId])
        {
            ++B1sparseId;
            ++mismatch;
        }
        else if(B1rowId[B1sparseId] > B2rowId[B2sparseId])
        {
            ++B2sparseId;
            ++mismatch;
        }
        else //Match
        {
            ++B1sparseId;
            ++B2sparseId;
        }
     }
  }
  return mismatch;
}


int sparseUpdateAndCheckPattern(const double tol, const int nsize, const int nnzA,
		   const double A[], const int ArowId[],
 	           const int AcolSum[], int *nnzB, double B[], 
		   int BrowId[], int BcolSum[], const double gamma,
                   bool doPatternCheck, bool strictSamePattern)
{
  int j,k,mismatch;

  std::vector<int> origBcolSum(nsize+1);
  std::vector<int> origBrowId(*nnzB);

  if(doPatternCheck)
  {
     // Copy the pattern of the old reduced matrix  
     for(j = 0; j<nsize; ++j) //iterating over columns
     {
        origBcolSum[j] = BcolSum[j];
        for(k = BcolSum[j]; k < BcolSum[j+1]; ++k) //iterating over rows
        {
            origBrowId[k] = BrowId[k];
        }
     }
     origBcolSum[nsize] = BcolSum[nsize];
  }

  //Re-compute thresholded matrix
  if(doPatternCheck && strictSamePattern)
  {
      sparseUpdateGamma(nsize, nnzA, A, ArowId, AcolSum,
                               *nnzB, B, BrowId, BcolSum, gamma);
      mismatch = 0;
      return mismatch;
  }
  else
  {
      sparseOffDiagThreshCopyGamma(tol, nsize,
	    nnzA, A, ArowId, AcolSum,
	    nnzB, B, BrowId, BcolSum, gamma);
  }

  mismatch = nsize*nsize;

  if(doPatternCheck)
  {
      mismatch = checkPattern(nsize, origBcolSum[nsize], &origBrowId[0], &origBcolSum[0],
                                     BcolSum[nsize], BrowId, BcolSum);
  }
  return mismatch;
}


#if defined SUNDIALS2
static int sparseToDense(void *user_data, DlsMat Jac)
{
  Jsparse *sm= ((cv_param *) user_data)->sparseMtx;

  for(int j=0; j < sm->nSize; ++j)
  {
    realtype* col = DENSE_COL(Jac,j);
    for(int k = sm->mtxColSum[j]; k < sm->mtxColSum[j+1]; ++k)
    {
      col[sm->mtxRowIdx[k]] = sm->mtxData[k];
    }
  }
  return 0;
}
#elif defined SUNDIALS3 || defined SUNDIALS4
static int sparseToDense(void *user_data, SUNMatrix Jac)
{
  Jsparse *sm= ((cv_param *) user_data)->sparseMtx;

  for(int j=0; j < sm->nSize; ++j)
  {
    for(int k = sm->mtxColSum[j]; k < sm->mtxColSum[j+1]; ++k)
    {
      int row = sm->mtxRowIdx[k];
      SM_ELEMENT_D(Jac,row,j) = sm->mtxData[k];
    }
  }
  return 0;
}
#endif

#if defined SUNDIALS2
int jac_full_dense(long int N, realtype t, N_Vector y, N_Vector fy,
			DlsMat Jac, void *user_data,
			N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
#elif defined SUNDIALS3 || defined SUNDIALS4
int jac_full_dense(realtype t, N_Vector y, N_Vector fy,
                   SUNMatrix Jac, void *user_data,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
#endif
{
  cv_param *cvp=(cv_param *)user_data;
  double startTime;

  int j,flag;
  if(cvp->sparseMtx->abort_matrix_eval) { return -1; }

  setupJacobianSparse(t,y,fy,cvp,tmp1,tmp2,tmp3);
  // At this point the entire Jacobian matrix should be constructed
  // in sparse form
  //
  
  //Setup Jacobian Sparse times itself 
  startTime = getHighResolutionTime();

  flag = sparseToDense(user_data,Jac);

  cvp->jacSetupTime += getHighResolutionTime() - startTime;
  ++(cvp->nJacSetup);

  return flag;
}
