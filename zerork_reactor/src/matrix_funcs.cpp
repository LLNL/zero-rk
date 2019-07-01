#include "matrix_funcs.h"
#include "cv_param_sparse.h"
#include "ode_funcs.h"
#include "utility_funcs.h"
#include "omp.h"

#include <limits>

static int get_fillfactor(int n, int* colSum, int* rowIdx, double* data, int * perm_c, superlu_options_t options);


int jac_full_prec_setup(realtype t, N_Vector y, N_Vector fy,
			booleantype jok, booleantype *jcurPtr,
			realtype gamma, void *user_data,
			N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
//  printf("ZERORK jac_full_prec_setup\n");
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
//  long int currNumErrTestFails;
//  CVodeGetNumErrTestFails(cvp->cvodeMemPtr, &currNumErrTestFails); 
//  CVodeGetNumNonlinSolvConvFails(cvp->cvodeMemPtr, &currNumErrTestFails);  
  if(cvp->sparseMtx->abort_matrix_eval) { return -1; }

  if(cvp->sparseMtx->isFirstFactor[0])
   {
      lunnz_last_perm = NV_LENGTH_S(y)*NV_LENGTH_S(y);
      ff_last_perm = 100.0;
      gamma_last_perm = 1.0;
   }


  //circos plots
#ifdef ZERORK_CIRCOS_PLOTS
  {
     char filename_prefix[64];
     snprintf(filename_prefix,64,"jac_out%03d",num_calls);
     double * y_ptr = NV_DATA_S(y);
     double Temp = y_ptr[cvp->nSpc]*cvp->Tref;
     double Pres = cvp->Press[cvp->currReactor];
     cvp->mech->generateCircosFilesTPY(Temp, Pres,
                                       y_ptr,
                                       30, //num_reacs
                                       1, //exclude_unimolecular
                                       filename_prefix);
  }
#endif


  if(!jok) // if the Jacobian is not ok, process a new jacobian
    {
      (*jcurPtr)=TRUE; //indicate that Jacobian data was recomputed
      setupJacobianSparse(t,y,fy,user_data,tmp1,tmp2,tmp3);
      // At this point the entire Jacobian matrix should be constructed
    }
  else //jok
    {
      (*jcurPtr)=FALSE; //indicate that Jacobian data was not recomputed
      if(cvp->sparseMtx->fakeUpdate) {return 0;}
    }

    startTime = getHighResolutionTime();
    bool bigOrderChange = fabs(log10(gamma/gamma_last_perm)) > cvp->sparseMtx->maxGammaChangeOrder;
    bool strictSamePattern = cvp->sparseMtx->strictSamePattern
                                      && jok
                                      && !bigOrderChange;
    bool doPatternCheck = (!(cvp->sparseMtx->isFirstFactor[0]) && cvp->sparseMtx->permThresh > 0.0);
    misMatch = sparseUpdateAndCheckPattern(cvp->sparseMtx->offDiagThreshold,
			    cvp->sparseMtx->nSize,
			    cvp->sparseMtx->nNonZero,
			    cvp->sparseMtx->mtxData,
			    cvp->sparseMtx->mtxRowIdx,
			    cvp->sparseMtx->mtxColSum,
			    &(cvp->sparseMtx->reduceNNZ[0]),
			    cvp->sparseMtx->reduceData,
			    cvp->sparseMtx->reduceRowIdx,
			    cvp->sparseMtx->reduceColSum,
                            gamma,cvp->sparseMtx->threshType,
                            doPatternCheck,
                            strictSamePattern);

    cvp->sparseMtx->optionSLU[0].Fact=DOFACT;
    if(
        strictSamePattern 
        || (
              //misMatch < (int)((double)cvp->sparseMtx->reduceNNZ * cvp->sparseMtx->permThresh)
              misMatch < (int)((double)cvp->sparseMtx->nSize * cvp->sparseMtx->permThresh)
              //cvp->sparseMtx->fillFactor < 1.01*ff_last_perm
              && cvp->sparseMtx->numPermReUses[0] < 40
              && !(cvp->sparseMtx->isFirstFactor[0])
              && !(bigOrderChange)
           )
      )
    {
         //if(cvp->sparseMtx->ILU) cvp->sparseMtx->optionSLU[0].Fact=SamePattern;
         cvp->sparseMtx->optionSLU[0].Fact=SamePattern;
         if(!strictSamePattern && jok) {++(cvp->sparseMtx->numPermReUses[0]);}
         ++(cvp->sparseMtx->numPermReUses[0]);
    } else {
//        printf("Doing col perm: jok = %d, rNNZ = %d, misMatch = %d, numReUses = %d, fillFactor = %g\n",jok,cvp->sparseMtx->reduceNNZ[0],misMatch,cvp->sparseMtx->numPermReUses[0], cvp->sparseMtx->fillFactor[0]);
        cvp->sparseMtx->numPermReUses[0] = 0;
    }

    cvp->precSetupTime += getHighResolutionTime() - startTime;

//  BEGIN: Debugging output
//          printf("preconditioner at %14.8e s\n",t);
//          print_sp_matrix(
//                          cvp->sparseMtx->nSize,cvp->sparseMtx->nSize,
//                          cvp->sparseMtx->reduceColSum,
//                          cvp->sparseMtx->reduceRowIdx,
//                          cvp->sparseMtx->reduceData
//                         );
//          printf("\n\n\n");
//  END:  Debugging output

//    printf("misMatch, jok, numReUses, netfails = %d, %d, %d, %d\n", misMatch, jok, numReUses, currNumErrTestFails);

    if(!(cvp->sparseMtx->isFirstFactor[0])  && cvp->sparseMtx->optionSLU[0].Fact != SamePattern_SameRowPerm)
    {
            Destroy_SuperNode_Matrix(&(cvp->sparseMtx->Lslu[0]));
            Destroy_CompCol_Matrix(&(cvp->sparseMtx->Uslu[0]));
    }

    //Update number of non-zeros in M.
    ((NCformat *)cvp->sparseMtx->Mslu[0].Store)->nnz = cvp->sparseMtx->reduceNNZ[0];

    cvp->sparseMtx->Bslu[0].ncol=0; // in dlinsolx1.c example this is supposed to
                                 // indicate that a solution is not needed only
                                 // the factorization

    // Perform permutation
    if(cvp->sparseMtx->optionSLU[0].Fact == DOFACT)
      {
        startTime = getHighResolutionTime();
        get_perm_c(MMD_AT_PLUS_A, &(cvp->sparseMtx->Mslu[0]), cvp->sparseMtx->colPermutation);
        ++(cvp->nColPerm); 
        cvp->colPermTime += getHighResolutionTime() - startTime;
      }


//    cvp->sparseMtx->optionSLU[0].PivotGrowth = YES;
    startTime = getHighResolutionTime();
    if(!cvp->sparseMtx->ILU)
      {
        dgssvx(&(cvp->sparseMtx->optionSLU[0]),&(cvp->sparseMtx->Mslu[0]),
	     cvp->sparseMtx->colPermutation,cvp->sparseMtx->rowPermutation,
	     cvp->sparseMtx->colElimTree,cvp->sparseMtx->equed,cvp->sparseMtx->Rvec,
	     cvp->sparseMtx->Cvec,&(cvp->sparseMtx->Lslu[0]),&(cvp->sparseMtx->Uslu[0]),
	     work,lwork,&(cvp->sparseMtx->Bslu[0]),&(cvp->sparseMtx->Xslu[0]),&rpg,&rcond,
	     ferr,berr,&(cvp->sparseMtx->mem_usage[0]),&(cvp->sparseMtx->statSLU[0]),
	     &flag);
      } else {
        dgsisx(&(cvp->sparseMtx->optionSLU[0]),&(cvp->sparseMtx->Mslu[0]),
	     cvp->sparseMtx->colPermutation,cvp->sparseMtx->rowPermutation,
	     cvp->sparseMtx->colElimTree,cvp->sparseMtx->equed,cvp->sparseMtx->Rvec,
	     cvp->sparseMtx->Cvec,&(cvp->sparseMtx->Lslu[0]),&(cvp->sparseMtx->Uslu[0]),
	     work,lwork,&(cvp->sparseMtx->Bslu[0]),&(cvp->sparseMtx->Xslu[0]),&rpg,&rcond,
	     &(cvp->sparseMtx->mem_usage[0]),&(cvp->sparseMtx->statSLU[0]),
	     &flag);
     }
    ++(cvp->nJacFactor);
//    cvp->sparseMtx->optionSLU[0].PivotGrowth = NO;
//    printf("Reciprocal pivot growth = %g\n",rpg);
//    printf("rpg, misMatch, jok = %g, %d, %d\n", rpg, misMatch, jok);

    cvp->sparseMtx->LUnnz[0] = ((SCformat *)cvp->sparseMtx->Lslu[0].Store)->nnz + ((NCformat *)cvp->sparseMtx->Uslu[0].Store)->nnz;
    cvp->sparseMtx->fillFactor[0] = (double)cvp->sparseMtx->LUnnz[0]/(double)cvp->sparseMtx->reduceNNZ[0];

    bool redoFact = false;

//    if(!(cvp->sparseMtx->isFirstFactor[0]) 
//        &&(
//                cvp->sparseMtx->fillFactor > 1.02*ff_last_perm
//             && cvp->sparseMtx->fillFactor > 2.0
//             && cvp->sparseMtx->optionSLU[0].Fact != DOFACT
//          )
//       )
//    {
//        redoFact = true;
//        numReUses = 0;
//    } 

//    if(rpg < 1.0e-8) redoFact = true;

    if(cvp->sparseMtx->optionSLU[0].DiagPivotThresh > 0.0)
    {
        for( j = 0; j < cvp->sparseMtx->nSize; ++j )
        {
            if( cvp->sparseMtx->rowPermutation[j] == -1 )
              {
                  int colElems = cvp->sparseMtx->reduceColSum[j+1] - cvp->sparseMtx->reduceColSum[j];
//                  printf("Negative value in row perm at %d.\nNelements in reduce matrix at column = %d\n\n",j,colElems);
                  redoFact = true;
              }
        }
    }

    if( redoFact ) 
    {
//         printf("Redoing factorization.  misMatch = %d, jok = %d\n", misMatch, jok);
              cvp->sparseMtx->optionSLU[0].Fact = DOFACT;
              Destroy_SuperNode_Matrix(&(cvp->sparseMtx->Lslu[0]));
              Destroy_CompCol_Matrix(&(cvp->sparseMtx->Uslu[0]));
              get_perm_c(MMD_AT_PLUS_A, &(cvp->sparseMtx->Mslu[0]), cvp->sparseMtx->colPermutation);
              dgssvx(&(cvp->sparseMtx->optionSLU[0]),&(cvp->sparseMtx->Mslu[0]),
	             cvp->sparseMtx->colPermutation,cvp->sparseMtx->rowPermutation,
	             cvp->sparseMtx->colElimTree,cvp->sparseMtx->equed,cvp->sparseMtx->Rvec,
	             cvp->sparseMtx->Cvec,&(cvp->sparseMtx->Lslu[0]),&(cvp->sparseMtx->Uslu[0]),
	             work,lwork,&(cvp->sparseMtx->Bslu[0]),&(cvp->sparseMtx->Xslu[0]),&rpg,&rcond,
	             ferr,berr,&(cvp->sparseMtx->mem_usage[0]),&(cvp->sparseMtx->statSLU[0]),
	             &flag);
        cvp->sparseMtx->LUnnz[0] = ((SCformat *)cvp->sparseMtx->Lslu[0].Store)->nnz + ((NCformat *)cvp->sparseMtx->Uslu[0].Store)->nnz;
        cvp->sparseMtx->fillFactor[0] = (double)cvp->sparseMtx->LUnnz[0]/(double)cvp->sparseMtx->reduceNNZ[0];
    }

    cvp->jacFactorTime += getHighResolutionTime() - startTime;

    if(cvp->sparseMtx->optionSLU[0].Fact == DOFACT)
    {
        ff_last_perm = cvp->sparseMtx->fillFactor[0];
        lunnz_last_perm = cvp->sparseMtx->LUnnz[0];
        gamma_last_perm = gamma;
    }

    if(cvp->sparseMtx->isFirstFactor[0]==1)
      {
        cvp->sparseMtx->isFirstFactor[0]=0;
      }

//   cvp->prevNumErrTestFails = currNumErrTestFails;
//    if(flag != 0)
//      {return flag;}
    // flag > 0, singular matrix, zero diagonal at row,col = flag
    // flag < 0, illegal input
    // positive return values 

  return flag;
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

  cvp->sparseMtx->optionSLU[0].Fact=FACTORED;
  cvp->sparseMtx->Bslu[0].ncol=1; // in dlinsolx1.c example this is reset to the number
                               // of right hand sides

  startTime = getHighResolutionTime();
  // copy the data to the sparse matrix RHS
  memcpy(cvp->sparseMtx->vecRHS,rptr,cpysize);

  //backsolve with expert driver function
  dgssvx(&(cvp->sparseMtx->optionSLU[0]),&(cvp->sparseMtx->Mslu[0]),
	 cvp->sparseMtx->colPermutation,cvp->sparseMtx->rowPermutation,
	 cvp->sparseMtx->colElimTree,cvp->sparseMtx->equed,cvp->sparseMtx->Rvec,
	 cvp->sparseMtx->Cvec,&(cvp->sparseMtx->Lslu[0]),&(cvp->sparseMtx->Uslu[0]),
	 work,lwork,&(cvp->sparseMtx->Bslu[0]),&(cvp->sparseMtx->Xslu[0]),&rpg,&rcond,
	 ferr,berr,&(cvp->sparseMtx->mem_usage[0]),&(cvp->sparseMtx->statSLU[0]),
	 &flag);


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
#if 1
       if(NV_Ith_S(tmp1,j) < cvp->minMassFrac)
         {NV_Ith_S(tmp1,j)=cvp->minMassFrac;}
#else
       if(fabs(tmp1_ptr[j]) < cvp->minMassFrac)
       {
         if(tmp1_ptr[j] < 0.0)
         {
             tmp1_ptr[j]= -cvp->minMassFrac;
         }
         else
         {
             tmp1_ptr[j]=  cvp->minMassFrac;
         }
       }
#endif
       NV_Ith_S(tmp2,j)=cvp->molWt[j]/(NV_Ith_S(tmp1,j))*cvp->invDens[0];
    }

    NV_Ith_S(tmp1,nSpc)=NV_Ith_S(y,nSpc);

    // calculate the reaction info for the strictly positive case
    startTimeDeriv = getHighResolutionTime();
    if(cvp->constPress)
    {
      const_dpdt_wsr(t,tmp1,tmp3,user_data);
    }
    else
    {
      const_vol_wsr(t,tmp1,tmp3,user_data);
    }
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
	sMptr[sparseId]+=(cvp->dTemp_dt[0])*(cvp->CvMass[j]*cvp->molWt[j]);
      }

    // step 3: divide by -1/(molWt[j]*meanCvMass)
    multFact=-1.0/(cvp->meanCvMass[0]);
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
    if(cvp->constPress)
    {
      const_dpdt_wsr(t,tmp1,tmp3,user_data);
    }
    else
    {
      const_vol_wsr(t,tmp1,tmp3,user_data);
    }
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
		   int BrowId[], int BcolSum[], const double gamma, const int threshType,
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
      if(threshType == 1)
        {
          sparseOffDiagThreshCopyGamma(tol, nsize,
			    nnzA, A, ArowId, AcolSum,
			    nnzB, B, BrowId, BcolSum,
                            gamma);
        } else if(threshType == 2) {
              sparseOffDiagThreshCopyGamma_cnorm(tol, nsize,
			    nnzA, A, ArowId, AcolSum,
			    nnzB, B, BrowId, BcolSum,
                            gamma);
        } else if(threshType == 3) {
              sparseOffDiagThreshCopyGamma_rnorm(tol, nsize,
			    nnzA, A, ArowId, AcolSum,
			    nnzB, B, BrowId, BcolSum,
                            gamma);
        } else if(threshType == 4) {
              sparseOffDiagThreshCopyGamma_rcmin(tol, nsize,
			    nnzA, A, ArowId, AcolSum,
			    nnzB, B, BrowId, BcolSum,
                            gamma);
        } else if(threshType == 5) {
              sparseOffDiagMagCopyGamma(tol, nsize,
			    nnzA, A, ArowId, AcolSum,
			    nnzB, B, BrowId, BcolSum,
                            gamma);
        } else {
           printf("Invalid threshType: %d.\n",threshType);
           exit(-1);
        }
  }

  mismatch = nsize*nsize;

  if(doPatternCheck)
  {
      mismatch = checkPattern(nsize, origBcolSum[nsize], &origBrowId[0], &origBcolSum[0],
                                     BcolSum[nsize], BrowId, BcolSum);
  }
  return mismatch;
}


//int sparse_jac_v(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, 	 
//                 void *user_data, N_Vector tmp)
//{
//  int i,j,ix;
//  cv_param *cvp=(cv_param *)user_data;
//  int nsize = cvp->sparseMtx->nSize;
//  double *sMptr=cvp->sparseMtx->mtxData;
//
//  //In brief testing, never called before PrecSetup
//  if(cvp->sparseMtx->isFirstFactor[0]==1)
//    {
//      return -1; 
//    }
//
//  //Re-compute jacobian
//  {
//    N_Vector tmp1, tmp2, tmp3;
//    tmp1 = N_VNew_Serial(nsize);
//    tmp2 = N_VNew_Serial(nsize);
//    tmp3 = N_VNew_Serial(nsize);
//    setupJacobianSparse(t,y,fy,cvp,tmp1,tmp2,tmp3);
//    N_VDestroy_Serial(tmp1);
//    N_VDestroy_Serial(tmp2);
//    N_VDestroy_Serial(tmp3);
//  }
//
//  N_VConst(0.0, Jv);
//
//  for (j = 0; j < nsize; ++j)
//    {
//      if (NV_Ith_S(v,j) != 0.)
//        {
//          for( i = cvp->sparseMtx->mtxColSum[j]; i < cvp->sparseMtx->mtxColSum[j+1]; ++i )
//            {
//              ix = cvp->sparseMtx->mtxRowIdx[i];
//              NV_Ith_S(Jv,ix) += NV_Ith_S(v,j) * sMptr[i];
//            }
//        }
//    }
//
//  return 0;
//}


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


int jac_full_dense(long int N, realtype t, N_Vector y, N_Vector fy,
			DlsMat Jac, void *user_data,
			N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  cv_param *cvp=(cv_param *)user_data;
  double startTime;

  int j,flag;
  if(cvp->sparseMtx->abort_matrix_eval) { return -1; }

  setupJacobianSparse(t,y,fy,user_data,tmp1,tmp2,tmp3);
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


//for interface to sparse direct solver
int jac_full_sparse(int N,
                    realtype t,
                    N_Vector y,
                    N_Vector fy,
                    double user_jac[],
                    void *user_data,
                    N_Vector tmp1,
                    N_Vector tmp2,
                    N_Vector tmp3)
{
  cv_param *cvp=(cv_param *)user_data;
  double startTime;

  int j,flag;
  if(cvp->sparseMtx->abort_matrix_eval) { return -1; }

  setupJacobianSparse(t,y,fy,user_data,tmp1,tmp2,tmp3);
  // At this point the entire Jacobian matrix should be constructed
  // in sparse form

  //Setup Jacobian Sparse times itself
  startTime = getHighResolutionTime();

  for(j=0; j < cvp->sparseMtx->nNonZero; ++j)
  {
    user_jac[j] = cvp->sparseMtx->mtxData[j];
  }

  cvp->jacSetupTime += getHighResolutionTime() - startTime;
  ++(cvp->nJacSetup);

  return 0;
}

// create a prototype for the internal cvode function not included in the
// public headers, but found in cvode_direct_impl.h
extern "C" int cvDlsDenseDQJac(long int N, realtype t,
                   N_Vector y, N_Vector fy,
                   DlsMat Jac, void *data,
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int jac_full_sparse_divided_diff(int N,
                    realtype t,
                    N_Vector y,
                    N_Vector fy,
                    double user_jac[],
                    void *user_data,
                    N_Vector tmp1,
                    N_Vector tmp2,
                    N_Vector tmp3)
{
  cv_param *cvp=(cv_param *)user_data;
  int flag, sparse_id;
  DlsMat diff_quotient=NewDenseMat(N,N);
  // compute the finite difference Jacobian
  flag = cvDlsDenseDQJac(N,
                         t,
                         y,
                         fy,
                         diff_quotient,
                         cvp->cvodeMemPtr,
                         tmp1,
                         tmp2,
                         tmp3);

  // dense to sparse
  sparse_id=0;
  for(int j=0; j<N; ++j) {
    realtype* dense_column = DENSE_COL(diff_quotient,j);
    for(int k=0; k<N; ++k) {
      if(k == cvp->sparseMtx->mtxRowIdx[sparse_id]) {
        user_jac[sparse_id] = dense_column[k];
        ++sparse_id;
      }
    }
  }

  // destroy dense matrix
  DestroyMat(diff_quotient);
  return 0;
}

