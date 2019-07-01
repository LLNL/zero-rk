#include "matrix_funcs.h"
#include "cv_param_sparse.h"
#include "ode_funcs.h"
#include "utility_funcs.h"

int get_fillfactor(int n, int* colSum, int* rowIdx, double* data, int * perm_c, superlu_options_t options);

int jac_full_prec_setup(realtype t, N_Vector y, N_Vector fy,
			booleantype jok, booleantype *jcurPtr,
			realtype gamma, void *user_data,
			N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  cv_param *cvp=(cv_param *)user_data;
  double startTime;

  static int numReUses = 0;
  int j;
  int flag;
  
  if(!jok) // if the Jacobian is not ok, process a new jacobian
    {
      setupJacobianSparse(t,y,fy,cvp,tmp1,tmp2,tmp3);
      
 /*     
      bool reThresh = true;
      if(!cvp->sparseMtx->isFirstFactor)
        {
          reThresh = reThreshCheck(cvp->sparseMtx->offDiagThreshold,
			    cvp->sparseMtx->nSize,
			    cvp->sparseMtx->nNonZero,
			    cvp->sparseMtx->mtxData,
			    cvp->sparseMtx->mtxRowIdx,
			    cvp->sparseMtx->mtxColSum,
			    cvp->sparseMtx->reduceNNZ,
			    cvp->sparseMtx->reduceData,
			    cvp->sparseMtx->reduceRowIdx,
			    cvp->sparseMtx->reduceColSum,
                            gamma);
        }
*/
      int misMatch=0;
      if(!cvp->sparseMtx->isFirstFactor)
        {
          misMatch = calcMismatch(cvp->sparseMtx->offDiagThreshold,
			    cvp->sparseMtx->nSize,
			    cvp->sparseMtx->nNonZero,
			    cvp->sparseMtx->mtxData,
			    cvp->sparseMtx->mtxRowIdx,
			    cvp->sparseMtx->mtxColSum,
			    cvp->sparseMtx->reduceNNZ,
			    cvp->sparseMtx->reduceData,
			    cvp->sparseMtx->reduceRowIdx,
			    cvp->sparseMtx->reduceColSum,
                            gamma);
//          printf("mismatch = %d\n",misMatch);
        }

      long int currNumErrTestFails;
      CVodeGetNumErrTestFails(cvp->cvodeMemPtr, &currNumErrTestFails); 
//        CVodeGetNumNonlinSolvConvFails(cvp->cvodeMemPtr, &currNumErrTestFails); 
//        if(cvp->sparseMtx->isFirstFactor == 1 || numReUses > 3)
//        if(reThresh || cvp->prevNumErrTestFails != currNumErrTestFails || numReUses > 10)
      if(1)
//      if(cvp->sparseMtx->isFirstFactor == 1 || cvp->prevNumErrTestFails != currNumErrTestFails || numReUses > 5)
        {
          numReUses = 0;
          cvp->prevNumErrTestFails = currNumErrTestFails;

          // copy only the off diagonal terms larger than the threshold
          if(cvp->sparseMtx->threshType == 1)
            {
              sparseOffDiagThreshCopyGamma(cvp->sparseMtx->offDiagThreshold,
  			    cvp->sparseMtx->nSize,
			    cvp->sparseMtx->nNonZero,
			    cvp->sparseMtx->mtxData,
			    cvp->sparseMtx->mtxRowIdx,
			    cvp->sparseMtx->mtxColSum,
			    &(cvp->sparseMtx->reduceNNZ),
			    cvp->sparseMtx->reduceData,
			    cvp->sparseMtx->reduceRowIdx,
			    cvp->sparseMtx->reduceColSum,
                            gamma);
            } else if(cvp->sparseMtx->threshType == 2) {
              sparseOffDiagThreshCopyGamma_cnorm(cvp->sparseMtx->offDiagThreshold,
			    cvp->sparseMtx->nSize,
			    cvp->sparseMtx->nNonZero,
			    cvp->sparseMtx->mtxData,
			    cvp->sparseMtx->mtxRowIdx,
			    cvp->sparseMtx->mtxColSum,
			    &(cvp->sparseMtx->reduceNNZ),
			    cvp->sparseMtx->reduceData,
			    cvp->sparseMtx->reduceRowIdx,
			    cvp->sparseMtx->reduceColSum,
                            gamma);
            } else if(cvp->sparseMtx->threshType == 3) {
              sparseOffDiagThreshCopyGamma_rnorm(cvp->sparseMtx->offDiagThreshold,
			    cvp->sparseMtx->nSize,
			    cvp->sparseMtx->nNonZero,
			    cvp->sparseMtx->mtxData,
			    cvp->sparseMtx->mtxRowIdx,
			    cvp->sparseMtx->mtxColSum,
			    &(cvp->sparseMtx->reduceNNZ),
			    cvp->sparseMtx->reduceData,
			    cvp->sparseMtx->reduceRowIdx,
			    cvp->sparseMtx->reduceColSum,
                            gamma);
            } else if(cvp->sparseMtx->threshType == 4) {
              sparseOffDiagThreshCopyGamma_rcmin(cvp->sparseMtx->offDiagThreshold,
			    cvp->sparseMtx->nSize,
			    cvp->sparseMtx->nNonZero,
			    cvp->sparseMtx->mtxData,
			    cvp->sparseMtx->mtxRowIdx,
			    cvp->sparseMtx->mtxColSum,
			    &(cvp->sparseMtx->reduceNNZ),
			    cvp->sparseMtx->reduceData,
			    cvp->sparseMtx->reduceRowIdx,
			    cvp->sparseMtx->reduceColSum,
                            gamma);
            } else if(cvp->sparseMtx->threshType == 5) {
              sparseOffDiagMagCopyGamma(cvp->sparseMtx->offDiagThreshold,
			    cvp->sparseMtx->nSize,
			    cvp->sparseMtx->nNonZero,
			    cvp->sparseMtx->mtxData,
			    cvp->sparseMtx->mtxRowIdx,
			    cvp->sparseMtx->mtxColSum,
			    &(cvp->sparseMtx->reduceNNZ),
			    cvp->sparseMtx->reduceData,
			    cvp->sparseMtx->reduceRowIdx,
			    cvp->sparseMtx->reduceColSum,
                            gamma);
            } else {
              printf("Invalid threshType: %d.\n",cvp->sparseMtx->threshType);
              return -1;
            }
          cvp->sparseMtx->optionSLU.Fact=DOFACT;
          if(!cvp->sparseMtx->isFirstFactor)
            {
    	      // You must destroy the L and U before each new factorization, or
	      // you will have a memory leak
	      Destroy_SuperNode_Matrix(&(cvp->sparseMtx->Lslu));
  	      Destroy_CompCol_Matrix(&(cvp->sparseMtx->Uslu));
              //if(cvp->sparseMtx->ILU) cvp->sparseMtx->optionSLU.Fact=SamePattern;
              if(misMatch < cvp->sparseMtx->nSize/3) cvp->sparseMtx->optionSLU.Fact=SamePattern; // 2ma opt: 4000
            }
        } else {
          ++numReUses;
          sparseUpdateGamma(cvp->sparseMtx->nSize,
			cvp->sparseMtx->nNonZero,
			cvp->sparseMtx->mtxData,
			cvp->sparseMtx->mtxRowIdx,
			cvp->sparseMtx->mtxColSum,
			cvp->sparseMtx->reduceNNZ,
			cvp->sparseMtx->reduceData,
			cvp->sparseMtx->reduceRowIdx,
			cvp->sparseMtx->reduceColSum,
                        gamma);
	    cvp->sparseMtx->optionSLU.Fact=SamePattern_SameRowPerm;
        }

      //Debugging printing output
      if(0)
        {
          //void print_sp_matrix(int m, int n, int* colSum, int* rowIdx, double* avals)
          printf("preconditioner at %14.8e s\n",t);
          print_sp_matrix(cvp->sparseMtx->nSize,cvp->sparseMtx->nSize,cvp->sparseMtx->reduceColSum,cvp->sparseMtx->reduceRowIdx,cvp->sparseMtx->reduceData);
          printf("\n\n\n");
        }
      (*jcurPtr)=TRUE; //indicate that Jacobian data was recomputed
    }
  else //jok
    {
      (*jcurPtr)=FALSE; //indicate that Jacobian data was not recomputed
      if(cvp->sparseMtx->fakeUpdate) {return 0;}

      sparseUpdateGamma(cvp->sparseMtx->nSize,
			cvp->sparseMtx->nNonZero,
			cvp->sparseMtx->mtxData,
			cvp->sparseMtx->mtxRowIdx,
			cvp->sparseMtx->mtxColSum,
			cvp->sparseMtx->reduceNNZ,
			cvp->sparseMtx->reduceData,
			cvp->sparseMtx->reduceRowIdx,
			cvp->sparseMtx->reduceColSum,
                        gamma);
      cvp->sparseMtx->optionSLU.Fact=SamePattern_SameRowPerm;
    }

    int lwork=0;
    void *work=NULL;
    double rpg,rcond;
    double ferr[1],berr[1]; // length is the number of RHS

    //dCreate_CompCol_Matrix(&(cvp->sparseMtx->Mslu),
    //			     cvp->sparseMtx->nSize,
    //			     cvp->sparseMtx->nSize,
    //			     cvp->sparseMtx->reduceNNZ,
    //			     cvp->sparseMtx->reduceData,
    //			     cvp->sparseMtx->reduceRowIdx,
    //			     cvp->sparseMtx->reduceColSum,
    //			     SLU_NC,SLU_D,SLU_GE);

    ((NCformat *)cvp->sparseMtx->Mslu.Store)->nnz=cvp->sparseMtx->reduceNNZ;

    cvp->sparseMtx->Bslu.ncol=0; // in dlinsolx1.c example this is supposed to
                                 // indicate that a solution is not needed only
                                 // the factorization

    // Perform permutation
    startTime = getHighResolutionTime();
    if(cvp->sparseMtx->optionSLU.Fact == DOFACT)
      {
         get_perm_c(MMD_AT_PLUS_A, &(cvp->sparseMtx->Mslu), cvp->sparseMtx->colPermutation);
      }
    ++(cvp->nColPerm); 
    cvp->colPermTime += getHighResolutionTime() - startTime;


    if(0)//if(t>1.0e-2)
      {
        printf("\n\n");
        printf("colPerm MMD_AT_PLUS_A : \n\n");
        for(j = 0; j < cvp->sparseMtx->nSize; ++j)
          {
             printf("%d\n",cvp->sparseMtx->colPermutation[j]);
          }
        printf("\n\n");
        printf("int xa[%d] = {\n",(cvp->sparseMtx->nSize+1));
        for(j = 0; j < cvp->sparseMtx->nSize; ++j)
          {
             printf("%d,\n",cvp->sparseMtx->reduceColSum[j]);
          }
        printf("%d};\n",cvp->sparseMtx->reduceColSum[cvp->sparseMtx->nSize]);
        printf("\n\n");
        printf("int asub[%d] = {\n", cvp->sparseMtx->reduceNNZ);
        for(j = 0; j < cvp->sparseMtx->reduceNNZ-1; ++j)
          {
             printf("%d,\n",cvp->sparseMtx->reduceRowIdx[j]);
          }
        printf("%d};\n",cvp->sparseMtx->reduceRowIdx[cvp->sparseMtx->reduceNNZ-1]);
        printf("\n\n");
        printf("double a[%d] = {\n",cvp->sparseMtx->reduceNNZ);
        for(j = 0; j < cvp->sparseMtx->reduceNNZ-1; ++j)
          {
             printf("%e,\n",cvp->sparseMtx->reduceData[j]);
          }
        printf("%e};\n",cvp->sparseMtx->reduceData[cvp->sparseMtx->reduceNNZ-1]);
        printf("\n\n");
        exit(-1);
      }


//#define MC_MIN_FILL
#ifdef MC_MIN_FILL
    if(cvp->sparseMtx->optionSLU.Fact == DOFACT)
    {
        int nTrials = 1000;
        int bestFill = cvp->sparseMtx->nSize*cvp->sparseMtx->nSize; // Completely full matrix
        int swp_a,swp_b,tmp_a;
        int perm_c[cvp->sparseMtx->nSize];
        swp_a = swp_b = 0; //Avoid compiler warning

        for( k = 0; k < nTrials; ++k)
        {
            for(j = 0; j < cvp->sparseMtx->nSize; ++j)
            {
               perm_c[j] = cvp->sparseMtx->colPermutation[j];
            }
            if( k > 0 )
            {
                swp_a = rand() % cvp->sparseMtx->nSize;
                swp_b = rand() % cvp->sparseMtx->nSize;
                if(swp_a == swp_b) continue;
                tmp_a = cvp->sparseMtx->colPermutation[swp_a];
                perm_c[swp_a] = perm_c[swp_b];
                perm_c[swp_b] = tmp_a;
            }

//            startTime = getHighResolutionTime();
            cvp->sparseMtx->LUnnz = get_fillfactor(cvp->sparseMtx->nSize,
                                                   cvp->sparseMtx->reduceColSum,
                                                   cvp->sparseMtx->reduceRowIdx,
                                                   cvp->sparseMtx->reduceData,
                                                   perm_c,
                                                   cvp->sparseMtx->optionSLU);
//            cvp->jacFactorTime += getHighResolutionTime() - startTime;
//            ++(cvp->nJacFactor);

            if( cvp->sparseMtx->LUnnz < bestFill )
            {
                bestFill = cvp->sparseMtx->LUnnz;
                for(j = 0; j < cvp->sparseMtx->nSize; ++j)
                {
                   cvp->sparseMtx->colPermutation[j] = perm_c[j];
                }
            }
            if(k % (nTrials/10) == 0) printf("%d | ",bestFill);
        }// for(k = 0; ...)
        printf("==\n");
    }
#endif

    startTime = getHighResolutionTime();
    if(!cvp->sparseMtx->ILU)
      {
        dgssvx(&(cvp->sparseMtx->optionSLU),&(cvp->sparseMtx->Mslu),
	     cvp->sparseMtx->colPermutation,cvp->sparseMtx->rowPermutation,
	     cvp->sparseMtx->colElimTree,cvp->sparseMtx->equed,cvp->sparseMtx->Rvec,
	     cvp->sparseMtx->Cvec,&(cvp->sparseMtx->Lslu),&(cvp->sparseMtx->Uslu),
	     work,lwork,&(cvp->sparseMtx->Bslu),&(cvp->sparseMtx->Xslu),&rpg,&rcond,
	     ferr,berr,&(cvp->sparseMtx->mem_usage),&(cvp->sparseMtx->statSLU),
	     &flag);
      } else {
        dgsisx(&(cvp->sparseMtx->optionSLU),&(cvp->sparseMtx->Mslu),
	     cvp->sparseMtx->colPermutation,cvp->sparseMtx->rowPermutation,
	     cvp->sparseMtx->colElimTree,cvp->sparseMtx->equed,cvp->sparseMtx->Rvec,
	     cvp->sparseMtx->Cvec,&(cvp->sparseMtx->Lslu),&(cvp->sparseMtx->Uslu),
	     work,lwork,&(cvp->sparseMtx->Bslu),&(cvp->sparseMtx->Xslu),&rpg,&rcond,
	     &(cvp->sparseMtx->mem_usage),&(cvp->sparseMtx->statSLU),
	     &flag);
     }
    ++(cvp->nJacFactor);

    cvp->sparseMtx->LUnnz = ((SCformat *)cvp->sparseMtx->Lslu.Store)->nnz + ((NCformat *)cvp->sparseMtx->Uslu.Store)->nnz;
    cvp->sparseMtx->fillFactor = cvp->sparseMtx->LUnnz/cvp->sparseMtx->reduceNNZ;

    cvp->jacFactorTime += getHighResolutionTime() - startTime;

    if(cvp->sparseMtx->isFirstFactor==1)
      {
        cvp->sparseMtx->isFirstFactor=0;
      }

    if(flag != 0)
      {return flag;}
    // flag > 0, singular matrix, zero diagonal at row,col = flag
    // flag < 0, illegal input
    // positive return values 

  return 0;
}


// Solve Pz = r where P in this case is the full numerical Jacobian
//   M=I-gamma*J
int jac_full_prec_solveV3(realtype t, N_Vector y, N_Vector fy,
		      N_Vector r, N_Vector z, realtype gamma,
		      realtype delta, int lr, void *user_data,
		      N_Vector tmp)
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
	 ferr,berr,&(cvp->sparseMtx->mem_usage),&(cvp->sparseMtx->statSLU),
	 &flag);

  //book keeping
  ++(cvp->nBackSolve);
  cvp->backsolveTime += getHighResolutionTime() - startTime;

  //copy solution to output
  memcpy(zptr,cvp->sparseMtx->vecSoln,cpysize);

  return flag;
}

void setupJacobianSparse(realtype t,N_Vector y,N_Vector fy,cv_param* cvp,N_Vector tmp1,N_Vector tmp2,N_Vector tmp3)
{
  double dTemp,RuTemp,multFact,startTime;
  int nSpc;
  int rowId,sparseId,sparseId2,nElem;
  int j,k,N;
  double *sMptr=cvp->sparseMtx->mtxData;
  void* user_data = (void*)cvp;

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
       if(NV_Ith_S(tmp1,j) < cvp->minMassFrac)
         {NV_Ith_S(tmp1,j)=cvp->minMassFrac;}
       NV_Ith_S(tmp2,j)=cvp->molWt[j]/(NV_Ith_S(tmp1,j)*cvp->Dens);
    }

    NV_Ith_S(tmp1,nSpc)=NV_Ith_S(y,nSpc);

    // calculate the reaction info for the strictly positive case
    const_vol_wsr(t,tmp1,tmp3,user_data);

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
    const_vol_wsr(t,tmp1,tmp3,user_data);

    // step 3: approximate d(ydot[k])/dT with finite difference
    sparseId=cvp->sparseMtx->mtxColSum[nSpc];
    for(k=0; k<N; k++)
      {
        sMptr[sparseId]=(NV_Ith_S(tmp3,k)-NV_Ith_S(fy,k))*multFact;
        sparseId++;
      }

    cvp->jacSetupTime += getHighResolutionTime() - startTime;
    ++(cvp->nJacSetup);
}



void sparseOffDiagThreshCopy(const double tol, const int nsize, const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int *nnzB, double B[], 
			     int BrowId[], int BcolSum[])
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
	  if(ArowId[AsparseId]==j || fabs(A[AsparseId]) > tol) 
	    { // directly copy diagonal and any off-diagonal terms larger than
	      // tol
	      B[BsparseId]=A[AsparseId];
	      BrowId[BsparseId]=ArowId[AsparseId];
	      BcolSum[j+1]++; // increment the column sum
	      BsparseId++;    // increment B sparse index 
	    }
	  AsparseId++; // increment the A sparse index
	}
    }
  (*nnzB)=BcolSum[nsize];
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

void sparseOffDiagThreshCopy_cnorm(const double tol, const int nsize, const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int *nnzB, double B[], 
			     int BrowId[], int BcolSum[])
{
  int j,k,nelem;
  int AsparseId=0;
  int BsparseId=0;
  double *Adiag;
 
  Adiag=(double *)malloc(sizeof(double)*nsize);
  AsparseId=0;

  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      k=0;
      AsparseId=AcolSum[j];
      while(k<nelem)
	{
	  if(ArowId[AsparseId]==j) 
	    {Adiag[j]=tol*fabs(A[AsparseId]); k=nelem;}
	  ++k;
	  ++AsparseId;
	}
    }

  AsparseId=BsparseId=0;
  BcolSum[0]=0;
  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      BcolSum[j+1]=BcolSum[j]; 
      for(k=0; k<nelem; k++)
	{
	  if(ArowId[AsparseId]==j || fabs(A[AsparseId]) > Adiag[j]) 
	    { // directly copy diagonal and any off-diagonal terms larger than
	      // tol
	      B[BsparseId]=A[AsparseId];
	      BrowId[BsparseId]=ArowId[AsparseId];
	      BcolSum[j+1]++; // increment the column sum
	      BsparseId++;    // increment B sparse index 
	    }
	  AsparseId++; // increment the A sparse index
	}
    }
  (*nnzB)=BcolSum[nsize];
  free(Adiag);
}

void sparseOffDiagThreshCopyGamma_cnorm(const double tol, const int nsize, const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int *nnzB, double B[], 
			     int BrowId[], int BcolSum[], const double gamma)
{
  int j,k,nelem;
  int AsparseId=0;
  int BsparseId=0;
  double *Adiag;
 
  Adiag=(double *)malloc(sizeof(double)*nsize);
  AsparseId=0;

  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      k=0;
      AsparseId=AcolSum[j];
      while(k<nelem)
	{
	  if(ArowId[AsparseId]==j) 
	    {Adiag[j]=tol*fabs(1.0-gamma*A[AsparseId]); k=nelem;}
	  ++k;
	  ++AsparseId;
	}
    }

  AsparseId=BsparseId=0;
  BcolSum[0]=0;
  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      BcolSum[j+1]=BcolSum[j]; 
      for(k=0; k<nelem; k++)
	{
          bool diag = ArowId[AsparseId] == j;
          double currVal = -gamma*A[AsparseId];
          if (diag) currVal += 1.0;
	  if(diag || fabs(currVal) > Adiag[j]) 
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
  free(Adiag);
}


void sparseOffDiagThreshCopy_rnorm(const double tol, const int nsize, const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int *nnzB, double B[], 
			     int BrowId[], int BcolSum[])
{
  int j,k,nelem;
  int AsparseId=0;
  int BsparseId=0;
  double *Adiag;
 
  Adiag=(double *)malloc(sizeof(double)*nsize);
  AsparseId=0;

  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      k=0;
      AsparseId=AcolSum[j];
      while(k<nelem)
	{
	  if(ArowId[AsparseId]==j) 
	    {Adiag[j]=tol*fabs(A[AsparseId]); k=nelem;}
	  ++k;
	  ++AsparseId;
	}
    }

  AsparseId=BsparseId=0;
  BcolSum[0]=0;
  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      BcolSum[j+1]=BcolSum[j]; 
      for(k=0; k<nelem; k++)
	{
	  if(ArowId[AsparseId]==j || fabs(A[AsparseId]) > Adiag[ArowId[AsparseId]]) 
	    { // directly copy diagonal and any off-diagonal terms larger than
	      // tol
	      B[BsparseId]=A[AsparseId];
	      BrowId[BsparseId]=ArowId[AsparseId];
	      BcolSum[j+1]++; // increment the column sum
	      BsparseId++;    // increment B sparse index 
	    }
	  AsparseId++; // increment the A sparse index
	}
    }
  (*nnzB)=BcolSum[nsize];
  free(Adiag);
}


void sparseOffDiagThreshCopyGamma_rnorm(const double tol, const int nsize, const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int *nnzB, double B[], 
			     int BrowId[], int BcolSum[], const double gamma)
{
  int j,k,nelem;
  int AsparseId=0;
  int BsparseId=0;
  double *Adiag;
 
  Adiag=(double *)malloc(sizeof(double)*nsize);
  AsparseId=0;

  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      k=0;
      AsparseId=AcolSum[j];
      while(k<nelem)
	{
	  if(ArowId[AsparseId]==j) 
	    {Adiag[j]=tol*fabs(1.0-gamma*A[AsparseId]); k=nelem;}
	  ++k;
	  ++AsparseId;
	}
    }

  AsparseId=BsparseId=0;
  BcolSum[0]=0;
  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      BcolSum[j+1]=BcolSum[j]; 
      for(k=0; k<nelem; k++)
	{
          bool diag = ArowId[AsparseId] == j;
          double currVal = -gamma*A[AsparseId];
          if (diag) currVal += 1.0;
	  if(diag || fabs(currVal) > Adiag[ArowId[AsparseId]]) 
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
  free(Adiag);
}

void sparseOffDiagThreshCopy_rcmin(const double tol, const int nsize, const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int *nnzB, double B[], 
			     int BrowId[], int BcolSum[])
{
  int j,k,nelem;
  int AsparseId=0;
  int BsparseId=0;
  double *Adiag;
 
  Adiag=(double *)malloc(sizeof(double)*nsize);
  AsparseId=0;

  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      k=0;
      AsparseId=AcolSum[j];
      while(k<nelem)
	{
	  if(ArowId[AsparseId]==j) 
	    {Adiag[j]=tol*fabs(A[AsparseId]); k=nelem;}
	  ++k;
	  ++AsparseId;
	}
    }

  AsparseId=BsparseId=0;
  BcolSum[0]=0;
  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      BcolSum[j+1]=BcolSum[j]; 
      for(k=0; k<nelem; k++)
	{
	  if(ArowId[AsparseId]==j ||
	     fabs(A[AsparseId]) > Adiag[ArowId[AsparseId]] ||
	     fabs(A[AsparseId]) > Adiag[j]) 
	    { // directly copy diagonal and any off-diagonal terms larger than
	      // tol
	      B[BsparseId]=A[AsparseId];
	      BrowId[BsparseId]=ArowId[AsparseId];
	      BcolSum[j+1]++; // increment the column sum
	      BsparseId++;    // increment B sparse index 
	    }
	  AsparseId++; // increment the A sparse index
	}
    }
  (*nnzB)=BcolSum[nsize];
  free(Adiag);
}

void sparseOffDiagThreshCopyGamma_rcmin(const double tol, const int nsize, const int nnzA,
			     const double A[], const int ArowId[],
			     const int AcolSum[], int *nnzB, double B[], 
			     int BrowId[], int BcolSum[], const double gamma)
{
  int j,k,nelem;
  int AsparseId=0;
  int BsparseId=0;
  double *Adiag;
 
  Adiag=(double *)malloc(sizeof(double)*nsize);
  AsparseId=0;

  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      k=0;
      AsparseId=AcolSum[j];
      while(k<nelem)
	{
	  if(ArowId[AsparseId]==j) 
	    {Adiag[j]=tol*fabs(1.0-gamma*A[AsparseId]); k=nelem;}
	  ++k;
	  ++AsparseId;
	}
    }

  AsparseId=BsparseId=0;
  BcolSum[0]=0;
  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      BcolSum[j+1]=BcolSum[j]; 
      for(k=0; k<nelem; k++)
	{
          bool diag = ArowId[AsparseId] == j;
          double currVal = -gamma*A[AsparseId];
          if (diag) currVal += 1.0;
	  if(diag ||
	     fabs(currVal) > Adiag[ArowId[AsparseId]] ||
	     fabs(currVal) > Adiag[j]) 
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
  free(Adiag);
}


int compare_double(const void *a, const void *b)
{
  double aval = *((double *)a);
  double bval = *((double *)b);
  return ((aval > bval) ? -1 : 1);
}

void sparseOffDiagMagCopy(const double frac, const int nsize, const int nnzA,
			  const double A[], const int ArowId[],
			  const int AcolSum[], int *nnzB, double B[], 
			  int BrowId[], int BcolSum[])
{
  int j,k,nelem,nkeep;
  int AsparseId=0;
  int BsparseId=0;
  double thresh;
  double *Asort;

 
  Asort=(double *)malloc(sizeof(double)*nnzA);

  for(j=0; j<nnzA; j++)
    {Asort[j]=fabs(A[j]);}

  qsort(Asort,nnzA,sizeof(double),compare_double);
//  sort(Asort,Asort+nnzA);
  nkeep=(int)(frac*(double)nnzA+0.5);
  thresh=Asort[nkeep-1];

  //printf("# threshold value = %.18g\n", thresh);
  //exit(-1);

  AsparseId=BsparseId=0;
  BcolSum[0]=0;
  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      BcolSum[j+1]=BcolSum[j]; 
      for(k=0; k<nelem; k++)
	{
	  if(ArowId[AsparseId]==j || fabs(A[AsparseId]) > thresh)
	    { // directly copy diagonal and any off-diagonal terms larger than
	      // tol
	      B[BsparseId]=A[AsparseId];
	      BrowId[BsparseId]=ArowId[AsparseId];
	      BcolSum[j+1]++; // increment the column sum
	      BsparseId++;    // increment B sparse index 
	    }
	  AsparseId++; // increment the A sparse index
	}
    }
  (*nnzB)=BcolSum[nsize];
  free(Asort);
}


void sparseOffDiagMagCopyGamma(const double frac, const int nsize, const int nnzA,
			  const double A[], const int ArowId[],
			  const int AcolSum[], int *nnzB, double B[], 
			  int BrowId[], int BcolSum[], const double gamma)
{
  int j,k,nelem,nkeep;
  int AsparseId=0;
  int BsparseId=0;
  double thresh;
  double *Asort;

  Asort=(double *)malloc(sizeof(double)*nnzA);

  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      AsparseId=AcolSum[j];
      for(k=0; k<nelem; ++k)
	{
           bool diag = ArowId[AsparseId] == j;
           double currVal = -gamma*A[AsparseId];
           if (diag) currVal += 1.0;
	   Asort[AsparseId]=fabs(currVal);
	   ++AsparseId;
	}
    }

  qsort(Asort,nnzA,sizeof(double),compare_double);
//  sort(Asort,Asort+nnzA);
  nkeep=(int)(frac*(double)nnzA+0.5);
  thresh=Asort[nkeep-1];

  AsparseId=BsparseId=0;
  BcolSum[0]=0;
  for(j=0; j<nsize; j++) // column number
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      BcolSum[j+1]=BcolSum[j]; 
      for(k=0; k<nelem; k++)
	{
          bool diag = ArowId[AsparseId] == j;
          double currVal = -gamma*A[AsparseId];
          if (diag) currVal += 1.0;
	  if(diag || fabs(currVal) > thresh)
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
  free(Asort);
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

bool reThreshCheck(const double tol, const int nsize, const int nnzA,
		   const double A[], const int ArowId[],
 	           const int AcolSum[], const int nnzB, const double B[], 
		   const int BrowId[], const int BcolSum[], const double gamma)
{
  int k,brow,bcol,arow,BsparseId,AsparseId,nElemAcol,lastAmatch;
  double maxDroppedVal = 0.0;
  int droppedVals = 0;
  int retainedVals = 0;

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
          double currVal = fabs(gamma*A[AsparseId]);
          if(brow != arow)
            {
	      maxDroppedVal = max(currVal,maxDroppedVal);
              if(currVal >= tol) ++droppedVals;
            } else {
              if(currVal <= tol) ++retainedVals;
              lastAmatch = k;
              if(BsparseId < BcolSum[bcol+1]) break; //not last val in B[:,j]
            }
         }
     }

  printf("droppedVals: %d, maxDroppedVal: %8.4e, retainedVals: %d, gamma: %8.4e\n",
          droppedVals,     maxDroppedVal,        retainedVals,     gamma);
//  if(droppedVals > nnzA/3 || maxDroppedVal > 25.0*tol || retainedVals > 500)
//  if(droppedVals > nnzA/3)
    {
      return true;
    }
  return false;
}


int calcMismatch(const double tol, const int nsize, const int nnzA,
		   const double A[], const int ArowId[],
 	           const int AcolSum[], const int nnzB, const double B[], 
		   const int BrowId[], const int BcolSum[], const double gamma)
{
  int j,k,nelem,nelemOld,mismatch;
  int AsparseId=0;
  int BsparseId=0;
//  int newB[nnzA];
  int newBcolSum[nsize];
//  int newBrowId[nnzA];
  
  AsparseId=BsparseId=0;
  newBcolSum[0]=0;
  for(j=0; j<nsize; j++) //
    {
      nelem=AcolSum[j+1]-AcolSum[j];
      newBcolSum[j+1]=newBcolSum[j]; 
      for(k=0; k<nelem; k++)
	{
          bool diag = ArowId[AsparseId] == j;
          double currVal = -gamma*A[AsparseId];
          if (diag) currVal += 1.0;
	  if(diag || fabs(currVal) > tol) 
	    { // directly copy diagonal and any off-diagonal terms larger than
	      // tol
//	      newB[BsparseId]=currVal;
//	      newBrowId[BsparseId]=ArowId[AsparseId];
	      newBcolSum[j+1]++; // increment the column sum
	      BsparseId++;    // increment B sparse index 
	    }
	  AsparseId++; // increment the A sparse index
	}
    }

  mismatch = 0;
  //generated new structure for reduced matrix on this step
  for(j = 0; j<nsize; ++j)
  {
     nelemOld = BcolSum[j+1] - BcolSum[j];  
     nelem = newBcolSum[j+1] - newBcolSum[j];  
     mismatch += abs(nelem-nelemOld);
  }
  return mismatch;
}



int sparse_jac_v(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, 	 
                 void *user_data, N_Vector tmp)
{
  int i,j,ix;
  cv_param *cvp=(cv_param *)user_data;
  int nsize = cvp->sparseMtx->nSize;
  double *sMptr=cvp->sparseMtx->mtxData;

  //In brief testing, never called before PrecSetup
  if(cvp->sparseMtx->isFirstFactor==1)
    {
      return -1; 
    }

  //Re-compute jacobian
  {
    N_Vector tmp1, tmp2, tmp3;
    tmp1 = N_VNew_Serial(nsize);
    tmp2 = N_VNew_Serial(nsize);
    tmp3 = N_VNew_Serial(nsize);
    setupJacobianSparse(t,y,fy,cvp,tmp1,tmp2,tmp3);
    N_VDestroy_Serial(tmp1);
    N_VDestroy_Serial(tmp2);
    N_VDestroy_Serial(tmp3);
  }

  N_VConst(0.0, Jv);

  for (j = 0; j < nsize; ++j)
    {
      if (NV_Ith_S(v,j) != 0.)
        {
          for( i = cvp->sparseMtx->mtxColSum[j]; i < cvp->sparseMtx->mtxColSum[j+1]; ++i )
            {
              ix = cvp->sparseMtx->mtxRowIdx[i];
              NV_Ith_S(Jv,ix) += NV_Ith_S(v,j) * sMptr[i];
            }
        }
    }

  return 0;
}


int get_fillfactor(int n, int* colSum, int* rowIdx, double* data, int * perm_c,  superlu_options_t options)
{
    int m,nnz,info;
    m = n; //Assume square matrix
    nnz = colSum[n];

    SuperLUStat_t stat;
    SuperMatrix A;
    int      *perm_r; /* row permutations from partial pivoting */
    SuperMatrix L;      /* factor L */
    SuperMatrix U;      /* factor U */
//    SCformat *Lstore;
//    NCformat *Ustore;

    dCreate_CompCol_Matrix(&A, m, n, nnz, data, rowIdx, colSum, SLU_NC, SLU_D, SLU_GE);

    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);

    int * work = NULL;
    int lwork = 0;
    int relax;
    int panel_size;
    int * etree;
    etree = intMalloc(A.ncol);
    SuperMatrix AC;

    panel_size = sp_ienv(1);
    relax = sp_ienv(2);

    if(options.ColPerm != MY_PERMC)  get_perm_c(options.ColPerm, &A, perm_c);

    sp_preorder(&options, &A, perm_c, etree, &AC);
    dgstrf(&options, &AC, relax, panel_size, etree, work, lwork, perm_c, perm_r,  &L, &U, &stat, &info);

    int lu_nnz = ((SCformat *)L.Store)->nnz + ((NCformat *)U.Store)->nnz;

    Destroy_SuperMatrix_Store(&A);
    Destroy_CompCol_Permuted(&AC);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

    StatFree(&stat);
    free(perm_r);
    free(etree);

    return lu_nnz;
}
