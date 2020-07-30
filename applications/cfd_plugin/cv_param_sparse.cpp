#include "cv_param_sparse.h"


JsparseTermList * alloc_JsparseTermList(const int inpFD, const int inpFC)
{
  JsparseTermList *w;

  w=(JsparseTermList *)malloc(sizeof(JsparseTermList));
  if(w==NULL)
    {
      printf("ERROR: allocation 1 failed in alloc_JsparseTermList(...)\n");
      return NULL;
    }
  w->nFwdDestroy=inpFD;
  w->nFwdCreate =inpFC;

  w->fwdDestroy=(JsparseIdx *)malloc(sizeof(JsparseIdx)*inpFD);
  if(w->fwdDestroy==NULL)
    {
      free(w);
      printf("ERROR: allocation 2 failed in alloc_JsparseTermList(...)\n");
      return NULL;
    }

  w->fwdCreate=(JsparseIdx *)malloc(sizeof(JsparseIdx)*inpFC);
  if(w->fwdCreate==NULL)
    {
      free(w->fwdDestroy);
      free(w);
      printf("ERROR: allocation 3 failed in alloc_JsparseTermList(...)\n");
      return NULL;
    }

  return w;
}

void free_JsparseTermList(JsparseTermList *w)
{
  if(w != NULL)
    {
      free(w->fwdCreate);
      free(w->fwdDestroy);

      free(w);
    }
}

Jsparse * alloc_Jsparse(const zerork::mechanism &mechInp)
{
  int nSpc=mechInp.getNumSpecies();
  int nStep=mechInp.getNumSteps();

  int nReac,nProd;
  int j,k,m,rowIdx,colIdx,nzAddr;

  int nFD,nFC;

  int *isNonZero; // local dense matrix indicating the non-zero terms

  Jsparse *w;
  w=(Jsparse *)malloc(sizeof(Jsparse));
  if(w==NULL)
    {
      printf("ERROR: allocation 1 failed in alloc_Jsparse(...)\n");
      fflush(stdout);
      return NULL;
    }

  countJacobianTerms(mechInp,&nFD,&nFC);
  w->termList=alloc_JsparseTermList(nFD,nFC);
  if(w->termList==NULL)
    {
      free(w);
      printf("ERROR: allocation 2 failed in alloc_Jsparse(...)\n");
      fflush(stdout);
      return NULL;
    }

  w->nSize=nSpc+1; // matrix size

  // allocate temporary arrays
  isNonZero = (int *)malloc(sizeof(int)*(w->nSize)*(w->nSize));

  // allocate memory that is related to the matrix dimension
  w->mtxColSum  = (int *)malloc(sizeof(int)*((w->nSize+1)));
  w->diagIdx    = (int *)malloc(sizeof(int)*(w->nSize));
  w->lastRowIdx = (int *)malloc(sizeof(int)*(w->nSize));

  // initialize the dense nonzero flags
  for(j=0; j<nSpc; j++) // process the first nSpc columns
    {
      for(k=0; k<nSpc; k++) // process the first nSpc rows
	{isNonZero[j*(w->nSize)+k]=0;}

      isNonZero[j*(w->nSize)+j]=1;    // mark the diagonal
      isNonZero[j*(w->nSize)+nSpc]=1; // mark the last row
    }
  for(k=0; k<(w->nSize); k++) // mark nSize rows in the last column
    {isNonZero[nSpc*(w->nSize)+k]=1;}

  // re-parse the system filling in the Jacobian term data
  // Jacobian = d ydot(k)/ dy(j)
  nFD=nFC=0;

  for(j=0; j<nStep; j++)
    {
      nReac=mechInp.getOrderOfStep(j);
      nProd=mechInp.getNumProductsOfStep(j);
      for(k=0; k<nReac; k++)
	{
	  colIdx=mechInp.getSpecIdxOfStepReactant(j,k); // species being perturbed
	
	  // forward destruction
	  for(m=0; m<nReac; m++)
	    {
	      rowIdx=mechInp.getSpecIdxOfStepReactant(j,m); // species being destroyed
	      isNonZero[colIdx*(w->nSize)+rowIdx]=1; // mark location in dense

	      w->termList->fwdDestroy[nFD].concIdx=colIdx;
	      w->termList->fwdDestroy[nFD].rxnIdx=j;
	      w->termList->fwdDestroy[nFD].sparseIdx=rowIdx;
	      ++nFD;
	    }
	  // forward creation
	  for(m=0; m<nProd; m++)
	    {
	      rowIdx=mechInp.getSpecIdxOfStepProduct(j,m); // species being created
	      isNonZero[colIdx*(w->nSize)+rowIdx]=1; // mark location in dense

	      w->termList->fwdCreate[nFC].concIdx=colIdx;
	      w->termList->fwdCreate[nFC].rxnIdx=j;
	      w->termList->fwdCreate[nFC].sparseIdx=rowIdx;
	      ++nFC;
	    }
	}
    }

  // check the jacobian term count
  if(nFD != w->termList->nFwdDestroy)
    {
      printf("ERROR: in alloc_Jterm_param()\n");
      printf("       nFD = %d != %d = w->nFwdDestroy\n",nFD,
	     w->termList->nFwdDestroy);
      exit(-1);
    }
  if(nFC != w->termList->nFwdCreate)
    {
      printf("ERROR: in alloc_Jterm_param()\n");
      printf("       nFC = %d != %d = w->nFwdCreate\n",nFC,
	     w->termList->nFwdCreate);
      exit(-1);
    }

   // count the number of nonzero terms in the dense matrix
   w->nNonZero=1; // start the count at one so it can still serve as a flag
   w->mtxColSum[0]=0;
   for(j=0; j<(w->nSize); j++) // process column j
     {
       for(k=0; k<(w->nSize); k++)
	 {
	   if(isNonZero[j*(w->nSize)+k]==1)
	     {
	       isNonZero[j*(w->nSize)+k]=w->nNonZero;
	       w->nNonZero++;
	     }
	 }
       // after counting column j store the running total in column j+1
       // of the column sum
       w->mtxColSum[j+1]=w->nNonZero-1;
     }
   // now at each nonzero term, isNonZero is storing the (address+1) in the
   // actual compressed column storage data array

   w->nNonZero--; // decrement the count

   // allocate matrix data
   w->mtxData=(double *)malloc(sizeof(double)*w->nNonZero);
   w->mtxRowIdx=(int *)malloc(sizeof(int)*w->nNonZero);
   for(j=0; j<w->nNonZero; j++)
     {w->mtxData[j]=0.0;} // initialize to zero for the print function

   // scan the the isNonZero array to determine the row indexes
   // and the special data addresses
   for(j=0; j<(w->nSize); j++) // process column j
     {
       for(k=0; k<(w->nSize); k++)  // row k
	 {
	   nzAddr=isNonZero[j*(w->nSize)+k];
	   if(nzAddr>=1)
	     {w->mtxRowIdx[nzAddr-1]=k;}
	 }
       // record the diagonal address
       w->diagIdx[j]    =isNonZero[j*(w->nSize)+j]-1;
       w->lastRowIdx[j]=isNonZero[(j+1)*(w->nSize)-1]-1;
     }

   // use the isNonZero array as a lookup to store the proper compressed
   // column data storage
   for(j=0; j<nFD; j++)
     {
       rowIdx=w->termList->fwdDestroy[j].sparseIdx;
       colIdx=w->termList->fwdDestroy[j].concIdx;
       nzAddr=isNonZero[colIdx*(w->nSize)+rowIdx];
       if(nzAddr==0)
	 {
	   printf("ERROR: %d term in fwdDestroy points to a zero J term\n",
		  j);
	   printf("       at J(%d,%d)\n",rowIdx,colIdx);
	   exit(-1);
	 }
       w->termList->fwdDestroy[j].sparseIdx=nzAddr-1; // reset to sparse addr
     }
   for(j=0; j<nFC; j++)
     {
       rowIdx=w->termList->fwdCreate[j].sparseIdx;
       colIdx=w->termList->fwdCreate[j].concIdx;
       nzAddr=isNonZero[colIdx*(w->nSize)+rowIdx];
       if(nzAddr==0)
	 {
	   printf("ERROR: %d term in fwdCreate points to a zero J term\n",
		  j);
	   printf("       at J(%d,%d)\n",rowIdx,colIdx);
	   exit(-1);
	 }
       w->termList->fwdCreate[j].sparseIdx=nzAddr-1; // reset to sparse addr
     }


   // Now the JsparseTermList should contain a consistent set of indexes
   // for processing the Jacobian term by term.  Note that each element
   // in the list is independent, so the lists of JsparseIdx can be sorted
   // to try to improve cache locality

   // sparse solver specific allocations
   // SuperLU
   w->isFirstFactor=1;
   w->numPermReUses=0;
   /* Set the default input options:
	options.Fact = DOFACT;
        options.Equil = YES;
    	options.ColPerm = COLAMD;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
	options.DiagPivotThresh = 1.0;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES;
   */
   set_default_options(&(w->optionSLU));

//   w->optionSLU.DiagPivotThresh = DiagPivotThresh;
   w->optionSLU.ColPerm = MY_PERMC;
   w->optionSLU.RowPerm = LargeDiag;
   w->optionSLU.Equil=YES;

   w->vecRHS=(double *)malloc(sizeof(double)*(w->nSize));
   w->vecSoln=(double *)malloc(sizeof(double)*(w->nSize));
   w->rowPermutation=(int *)malloc(sizeof(int)*(w->nSize));
   w->colPermutation=(int *)malloc(sizeof(int)*(w->nSize));
   w->colElimTree=(int *)malloc(sizeof(int)*(w->nSize));

   w->fillFactor=1.0;
   w->LUnnz=w->nNonZero;
   w->reduceNNZ=w->nNonZero;
   w->reduceData=(double *)malloc(sizeof(double)*w->nNonZero);
   w->reduceRowIdx=(int *)malloc(sizeof(int)*w->nNonZero);
   w->reduceColSum=(int *)malloc(sizeof(int)*(w->nSize+1));

   // create a compressed column matrix
   // SLU_NC ==
   // SLU_D  == data type double precision
   // SLU_GE == matrix structure general
   //dCreate_CompCol_Matrix(&(w->Mslu),w->nSize,w->nSize,w->nNonZero,
   //			  w->mtxData,w->mtxRowIdx,w->mtxColSum,
   //			  SLU_NC,SLU_D,SLU_GE);
   w->Rvec=(double *)malloc(sizeof(double)*(w->nSize));
   w->Cvec=(double *)malloc(sizeof(double)*(w->nSize));

   dCreate_Dense_Matrix(&(w->Bslu),w->nSize,1,w->vecRHS,w->nSize, SLU_DN,SLU_D,SLU_GE);
   dCreate_Dense_Matrix(&(w->Xslu),w->nSize,1,w->vecSoln,w->nSize, SLU_DN,SLU_D,SLU_GE);

   dCreate_CompCol_Matrix(&(w->Mslu), w->nSize, w->nSize, w->nNonZero,
			  w->reduceData, w->reduceRowIdx,
			  w->reduceColSum, SLU_NC,SLU_D,SLU_GE);

   StatInit(&(w->statSLU)); // initialize SuperLU statistics
   w->abort_matrix_eval = false;

   free(isNonZero);
   return w;
}

void free_Jsparse(Jsparse *w)
{
  if(w!=NULL) {
    StatFree(&(w->statSLU));
    Destroy_SuperMatrix_Store(&(w->Bslu));
    Destroy_SuperMatrix_Store(&(w->Xslu));
    if(!w->isFirstFactor) {
      Destroy_SuperNode_Matrix(&(w->Lslu));
      Destroy_CompCol_Matrix(&(w->Uslu));
    }
    Destroy_SuperMatrix_Store(&(w->Mslu));


    free(w->vecRHS);
    free(w->vecSoln);

    free(w->reduceData);
    free(w->reduceRowIdx);
    free(w->reduceColSum);

    free(w->mtxData);
    free(w->mtxRowIdx);
    free(w->mtxColSum);
    free(w->diagIdx);
    free(w->lastRowIdx);

    // SuperLU data
    free(w->Rvec);
    free(w->Cvec);
    free(w->rowPermutation);
    free(w->colPermutation);
    free(w->colElimTree);

    free_JsparseTermList(w->termList);
    free(w);
  }
}

void countJacobianTerms(const zerork::mechanism &mechInp, int *nFwdDestroy, int *nFwdCreate)
{
  int j;
  int nStep=mechInp.getNumSteps();
  double nuSumReac,nuSumProd;

  (*nFwdDestroy)=0;
  (*nFwdCreate)=0;

  for(j=0; j<nStep; j++)
    {
      nuSumReac=mechInp.getOrderOfStep(j);
      nuSumProd=mechInp.getNumProductsOfStep(j);
      (*nFwdDestroy)+=nuSumReac*nuSumReac;
      (*nFwdCreate)+=nuSumReac*nuSumProd;
    }
}

double getElement_Jsparse(Jsparse *w,const int rowIdx, const int colIdx)
{
  int nElem=w->mtxColSum[colIdx+1]-w->mtxColSum[colIdx];
  int j;
  int currentRow;

  for(j=0; j<nElem; j++)
    {
      currentRow=w->mtxRowIdx[w->mtxColSum[colIdx]+j];
      if(rowIdx == currentRow)
	{return w->mtxData[w->mtxColSum[colIdx]+j];}
      else if(rowIdx < currentRow)
	{return 0.0;}
    }
  return 0.0;
}

void calcReaction_Jsparse(Jsparse *w, const double invPosConc[],
			  const double fwdROP[])
{
  int j;
  int rxnId,concId,sparseId;

  // set the full sparse array
  for(j=0; j<w->nNonZero; j++)
    {w->mtxData[j]=0.0;}

  // process the forward destruction terms
  for(j=0; j<w->termList->nFwdDestroy; j++)
    {
      concId   = w->termList->fwdDestroy[j].concIdx;
      rxnId    = w->termList->fwdDestroy[j].rxnIdx;
      sparseId = w->termList->fwdDestroy[j].sparseIdx;
      //printf(" j = %d, concId = %d, rxnId = %d, sparseId = %d\n",j,concId,rxnId,sparseId); fflush(stdout);
     w->mtxData[sparseId]-=fwdROP[rxnId]*invPosConc[concId];
    }
  // process the forward creation terms
  for(j=0; j<w->termList->nFwdCreate; j++)
    {
      concId   = w->termList->fwdCreate[j].concIdx;
      rxnId    = w->termList->fwdCreate[j].rxnIdx;
      sparseId = w->termList->fwdCreate[j].sparseIdx;
      w->mtxData[sparseId]+=fwdROP[rxnId]*invPosConc[concId];
    }
}

void change_JsparseThresh(Jsparse *w, double newThresh)
{w->offDiagThreshold=newThresh;}


