#include "cv_param_sparse.h"


static void transposeSparseMtxIndices(Jsparse *w);

//n is length of perm_c which is allocated by calling function
static void read_perm_c_from_file(int n, int* perm_c)
{
    int ncheck = 0;
    char line[1024];
    FILE *infile;
    infile=fopen("perm_c.txt", "r");

    //assumes perm_c has correct number of lines
    while (fgets(line, 1023, infile) != NULL && ncheck <= n) {
        perm_c[ncheck++] = atoi(line);
    }
    if(ncheck != n)
      {
        printf("ERROR: perm_c.txt not consistent with mechanism.");
        exit(-1);
      }
    return;
}


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

Jsparse * alloc_Jsparse(const zerork::mechanism &mechInp, int nReactors)
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
  w->nReactors=nReactors;

  w->nSize=nSpc+1; // matrix size

  // allocate temporary arrays
  isNonZero = (int *)malloc(sizeof(int)*(w->nSize)*(w->nSize));

  // allocate memory that is related to the matrix dimension
  w->mtxColSum  = (int *)malloc(sizeof(int)*((w->nSize*nReactors+1)));
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
//   printf("# number of nonzero terms = %d\n",w->nNonZero);
   fflush(stdout);

   // allocate matrix data
   w->mtxData=(double *)malloc(sizeof(double)*w->nNonZero*nReactors);
   w->mtxRowIdx=(int *)malloc(sizeof(int)*w->nNonZero*nReactors);
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
       w->termList->fwdDestroy[j].sparseIdxCSR=rowIdx; // save for CSR
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
       w->termList->fwdCreate[j].sparseIdxCSR=rowIdx; // save for CSR
       w->termList->fwdCreate[j].sparseIdx=nzAddr-1; // reset to sparse addr
     }


   // Now do the same for CSR
   w->mtxRowSum  = (int *)malloc(sizeof(int)*(w->nSize*nReactors+1));

   // count the number of nonzero terms in the dense matrix
   w->nNonZero=1; // start the count at one so it can still serve as a flag
   w->mtxRowSum[0]=0;
   for(j=0; j<(w->nSize); j++) // process row j
     {
       for(k=0; k<(w->nSize); k++)
	 {
	   if(isNonZero[k*(w->nSize)+j]>0)
	     {
	       isNonZero[k*(w->nSize)+j]=w->nNonZero;
	       w->nNonZero++;
	     }
	 }
       // after counting row j store the running total in row j+1
       // of the column sum
       w->mtxRowSum[j+1]=w->nNonZero-1;
     }

   w->nNonZero--; // decrement the count

   // allocate matrix data
   w->mtxColIdx=(int *)malloc(sizeof(int)*w->nNonZero*nReactors);

   // scan the the isNonZero array to determine the row indexes
   // and the special data addresses
   for(j=0; j<(w->nSize); j++) // process row j
     {
       for(k=0; k<(w->nSize); k++)  // col k
	 {
	   nzAddr=isNonZero[k*(w->nSize)+j];
	   if(nzAddr>0)
	     {w->mtxColIdx[nzAddr-1]=k;}
	 }
     }

   // use the isNonZero array as a lookup to store the proper compressed
   // column data storage
   for(j=0; j<nFD; j++)
     {
       rowIdx=w->termList->fwdDestroy[j].sparseIdxCSR;
       colIdx=w->termList->fwdDestroy[j].concIdx;
       nzAddr=isNonZero[colIdx*(w->nSize)+rowIdx];
       if(nzAddr==0)
	 {
	   printf("ERROR: %d term in fwdDestroy points to a zero J term\n",
		  j);
	   printf("       at J(%d,%d)\n",rowIdx,colIdx);
	   exit(-1);
	 }
       w->termList->fwdDestroy[j].sparseIdxCSR=nzAddr-1; // reset to sparse addr
     }
   for(j=0; j<nFC; j++)
     {
       rowIdx=w->termList->fwdCreate[j].sparseIdxCSR;
       colIdx=w->termList->fwdCreate[j].concIdx;
       nzAddr=isNonZero[colIdx*(w->nSize)+rowIdx];
       if(nzAddr==0)
	 {
	   printf("ERROR: %d term in fwdCreate points to a zero J term\n",
		  j);
	   printf("       at J(%d,%d)\n",rowIdx,colIdx);
	   exit(-1);
	 }
       w->termList->fwdCreate[j].sparseIdxCSR=nzAddr-1; // reset to sparse addr
     }


   // Now the JsparseTermList should contain a consistent set of indexes
   // for processing the Jacobian term by term.  Note that each element
   // in the list is independent, so the lists of JsparseIdx can be sorted
   // to try to improve cache locality

   // sparse solver specific allocations
   // SuperLU
   w->isFirstFactor=(int *)malloc(sizeof(int)*nReactors);
   for(k=0;k<nReactors;++k)
   {
       w->isFirstFactor[k]=1; // indicate to do the full factorization
   }
   w->numPermReUses=(int *)malloc(sizeof(int)*nReactors);
   for(k=0;k<nReactors;++k)
   {
       w->numPermReUses[k]=0; // indicate to do the full factorization
   }
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
   w->optionSLU=(superlu_options_t*)malloc(sizeof(superlu_options_t)*nReactors);
   for(k=0;k<nReactors;++k)
   {
       set_default_options(&(w->optionSLU[k]));

   /* Set the default options for ILU
        options.Fact = DOFACT;
        options.Equil = YES;
        options.ColPerm = COLAMD;
        options.Trans = NOTRANS;
        options.IterRefine = NOREFINE;
        options.DiagPivotThresh = 0.1; //different from complete LU
        options.SymmetricMode = NO;
        options.PivotGrowth = NO;
        options.ConditionNumber = NO;
        options.PrintStat = YES;
        options.RowPerm = LargeDiag;
        options.ILU_DropTol = 1e-4;
        options.ILU_FillTol = 1e-2;
        options.ILU_FillFactor = 10.0;
        options.ILU_DropRule = DROP_BASIC | DROP_AREA;
        options.ILU_Norm = INF_NORM;
        options.ILU_MILU = SILU;
   */
       ilu_set_default_options(&(w->optionSLU[k]));

       w->optionSLU[k].ILU_DropTol = 5e-4;
       w->optionSLU[k].ILU_FillTol = 1e-4;
       w->optionSLU[k].ILU_FillFactor = 2.0;
//       w->optionSLU[k].DiagPivotThresh = DiagPivotThresh;
       w->optionSLU[k].ILU_DropRule = DROP_BASIC | DROP_AREA;
       w->optionSLU[k].ILU_MILU = SILU; //SMILU_{1,2,3};
//       w->optionSLU[k].ColPerm = MMD_AT_PLUS_A; // best for iso-octane
//       w->optionSLU[k].ColPerm = COLAMD;
//       w->optionSLU[k].ColPerm = MMD_ATA;
       w->optionSLU[k].ColPerm = MY_PERMC;
       w->optionSLU[k].RowPerm = LargeDiag;
       w->optionSLU[k].Equil=YES;
//       w->optionSLU[k].SymmetricMode = YES; //Small DiagPivotThresh and MMD_AT_PLUS_A
//       w->optionSLU[k].PivotGrowth = YES;
   }

   w->equed=(char *)malloc(sizeof(char)*nReactors);
   w->vecRHS=(double *)malloc(sizeof(double)*(w->nSize*nReactors));
   w->vecSoln=(double *)malloc(sizeof(double)*(w->nSize*nReactors));
   w->rowPermutation=(int *)malloc(sizeof(int)*(w->nSize)*nReactors);
   w->colPermutation=(int *)malloc(sizeof(int)*(w->nSize)*nReactors);
   w->colElimTree=(int *)malloc(sizeof(int)*(w->nSize)*nReactors);

   w->fillFactor=(double *)malloc(sizeof(double)*nReactors);
   w->LUnnz=(int *)malloc(sizeof(int)*nReactors);
   w->reduceNNZ=(int *)malloc(sizeof(int)*nReactors);
   for(k=0;k<nReactors;++k)
   {
       w->reduceNNZ[k]=w->nNonZero;
   }
   w->reduceData=(double *)malloc(sizeof(double)*w->nNonZero*nReactors);
   w->reduceRowIdx=(int *)malloc(sizeof(int)*w->nNonZero*nReactors);
   w->reduceColSum=(int *)malloc(sizeof(int)*(w->nSize+1)*nReactors);
   w->reduceColIdx=(int *)malloc(sizeof(int)*w->nNonZero*nReactors);//CSR format
   w->reduceRowSum=(int *)malloc(sizeof(int)*(w->nSize+1)*nReactors); //CSR format
   w->transform_idx=(int *)malloc(sizeof(int)*w->nNonZero); //CSR format

   // create a compressed column matrix
   // SLU_NC ==
   // SLU_D  == data type double precision
   // SLU_GE == matrix structure general
   //dCreate_CompCol_Matrix(&(w->Mslu),w->nSize,w->nSize,w->nNonZero,
   //			  w->mtxData,w->mtxRowIdx,w->mtxColSum,
   //			  SLU_NC,SLU_D,SLU_GE);
   w->Rvec=(double *)malloc(sizeof(double)*(w->nSize)*nReactors);
   w->Cvec=(double *)malloc(sizeof(double)*(w->nSize)*nReactors);

   w->Bslu = (SuperMatrix *)malloc(sizeof(SuperMatrix)*nReactors);
   w->Xslu = (SuperMatrix *)malloc(sizeof(SuperMatrix)*nReactors);

   for(k=0;k<nReactors;++k)
   {
       dCreate_Dense_Matrix(&(w->Bslu[k]),w->nSize,1,&(w->vecRHS[(w->nSize)*k]),w->nSize,
			SLU_DN,SLU_D,SLU_GE);
       dCreate_Dense_Matrix(&(w->Xslu[k]),w->nSize,1,&(w->vecSoln[(w->nSize)*k]),w->nSize,
			SLU_DN,SLU_D,SLU_GE);
   }

   w->Mslu = (SuperMatrix *)malloc(sizeof(SuperMatrix)*nReactors);
   w->Lslu = (SuperMatrix *)malloc(sizeof(SuperMatrix)*nReactors);
   w->Uslu = (SuperMatrix *)malloc(sizeof(SuperMatrix)*nReactors);
   for(k=0;k<nReactors;++k)
   {
       dCreate_CompCol_Matrix(&(w->Mslu[k]), w->nSize, w->nSize, w->nNonZero,
			     &w->reduceData[k*w->nNonZero], &w->reduceRowIdx[k*w->nNonZero],
			     &w->reduceColSum[k*(w->nSize+1)], SLU_NC,SLU_D,SLU_GE);
   }


   w->statSLU=(SuperLUStat_t*)malloc(sizeof(SuperLUStat_t)*nReactors);
   for(k=0;k<nReactors;++k)
   {
     StatInit(&(w->statSLU[k])); // initialize SuperLU statistics
   }

   w->mem_usage = (mem_usage_t*)malloc(sizeof(mem_usage_t)*nReactors);

#ifdef CVIDT_USE_METIS
   METIS_SetDefaultOptions(w->metis_options);

   //Set to -1 for default
   w->metis_options[METIS_OPTION_NUMBERING] = 0; // C-numbering
   w->metis_options[METIS_OPTION_PFACTOR] = 0; // Default: 0 ["Good values often in the range of 60 to 200"]
   w->metis_options[METIS_OPTION_CCORDER] = 1;  // Default: 0
   w->metis_options[METIS_OPTION_COMPRESS] = 1; // Default: 1
   w->metis_options[METIS_OPTION_NITER] = 10; // Default: 10
   w->metis_options[METIS_OPTION_NSEPS] = 1; // Default: 1
   w->metis_options[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP1SIDED; // Default: ??
//METIS_RTYPE_FM FM-based cut refinnement. (SegFault)
//METIS_RTYPE_GREEDY Greedy-based cut and volume refinement. (SegFault)
//METIS_RTYPE_SEP2SIDED Two-sided node FM refinement.
//METIS_RTYPE_SEP1SIDED One-sided node FM refinement.
   w->metis_options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM; // Default: SHEM
//METIS_CTYPE_RM Random matching.
//METIS_CTYPE_SHEM Sorted heavy-edge matching.
   w->metis_options[METIS_OPTION_UFACTOR] = 1; // Default: 1 or 30 (ptype=rb or ptype=kway)
   w->metis_options[METIS_OPTION_SEED] = 0; // Default: ??
#endif

   w->reactorErrMult = (double*) malloc(sizeof(double)*nReactors);

   w->abort_matrix_eval = false;

   free(isNonZero);
   return w;
}

void free_Jsparse(Jsparse *w)
{
  if(w!=NULL)
    {
      for(int k=0; k<w->nReactors;++k)
      {
          StatFree(&(w->statSLU[k]));
          Destroy_SuperMatrix_Store(&(w->Bslu[k]));
          Destroy_SuperMatrix_Store(&(w->Xslu[k]));
          if(!w->isFirstFactor[k])
          {
              Destroy_SuperNode_Matrix(&(w->Lslu[k]));
              Destroy_CompCol_Matrix(&(w->Uslu[k]));
          }
//          Destroy_CompCol_Matrix(&(w->Mslu[k]));
          Destroy_SuperMatrix_Store(&(w->Mslu[k]));
      }


      free(w->vecRHS);
      free(w->vecSoln);

      free(w->reduceData);
      free(w->reduceRowIdx);
      free(w->reduceColSum);
      free(w->reduceColIdx);  //CSR format
      free(w->reduceRowSum); //CSR format
      free(w->transform_idx); //CSR format

      free(w->Bslu);
      free(w->Xslu);

      free(w->Mslu);
      free(w->Lslu);
      free(w->Uslu);

      free(w->equed);
      free(w->isFirstFactor);
      free(w->numPermReUses);
      free(w->fillFactor);
      free(w->LUnnz);
      free(w->reduceNNZ);

      free(w->mtxData);
      free(w->mtxRowIdx);
      free(w->mtxColSum);
      free(w->mtxColIdx); //CSR format
      free(w->mtxRowSum); //CSR format
      free(w->diagIdx);
      free(w->lastRowIdx);

      // SuperLU data
      free(w->Rvec);
      free(w->Cvec);
      free(w->rowPermutation);
      free(w->colPermutation);
      free(w->colElimTree);
      free(w->optionSLU);
      free(w->statSLU);
      free(w->mem_usage);

      free(w->reactorErrMult);

      free_JsparseTermList(w->termList);
      free(w);
    }
}




// void print_Jsparse(Jsparse *w)
// {
//   int j,k,nElems;
//   int rowIdx;
//   double val;
//   printf("# Number of nonzero elements: %d\n",w->nNonZero);

//   for(j=0; j<w->nSize; j++)
//     {
//       nElems=w->mtxColSum[j+1]-w->mtxColSum[j];
//       printf("   Col %d has %d nonzero elements\n",j,nElems);
//       for(k=0; k<nElems; k++)
// 	{
// 	  val=w->mtxData[w->mtxColSum[j]+k];
// 	  rowIdx=w->mtxRowIdx[w->mtxColSum[j]+k];
// 	  printf("     J(%d,%d) = %14.6e\n",rowIdx,j,val);
// 	}
//       printf("\n");
//       //exit(-1);
//     }

//   printf("# Number of forward destruction terms: %d\n",
// 	 w->termList->nFwdDestroy);
//   printf("# Number of forward creation terms: %d\n",
// 	 w->termList->nFwdCreate);
//   printf("# Number of reverse destruction terms: %d\n",
// 	 w->termList->nRevDestroy);
//   printf("# Number of reverse creation terms: %d\n",
// 	 w->termList->nRevCreate);
//   printf("# Forward Destruction Terms:\n");
//   printf("# i.e. both species are reactants in the same reaction\n");
//   printf("#\n#                  numer   denom\n");
//   printf("# J term   rxn ID  spc ID  spc ID  sparse mtx ID\n");
//   for(j=0; j<(w->termList->nFwdDestroy); j++)
//     {
//       printf("%6d  %6d  %6d  %6d  %6d\n",j,
// 	     w->termList->fwdDestroy[j].rxnIdx,
// 	     w->mtxRowIdx[w->termList->fwdDestroy[j].sparseIdx],
// 	     w->termList->fwdDestroy[j].concIdx,
// 	     w->termList->fwdDestroy[j].sparseIdx);
//     }
//   printf("\n\n");
//   printf("# Forward Creation Terms:\n");
//   printf("# i.e. numerator species is a product and the denominator\n");
//   printf("#      is a reactant in the same reaction\n");
//   printf("#\n#                  numer   denom\n");
//   printf("# J term   rxn ID  spc ID  spc ID\n");
//   for(j=0; j<(w->termList->nFwdCreate); j++)
//     {
//       printf("%6d  %6d  %6d  %6d  %6d\n",j,
// 	     w->termList->fwdCreate[j].rxnIdx,
// 	     w->mtxRowIdx[w->termList->fwdCreate[j].sparseIdx],
// 	     w->termList->fwdCreate[j].concIdx,
// 	     w->termList->fwdCreate[j].sparseIdx);
//     }
//   printf("\n\n");
//   printf("# Reverse Destruction Terms:\n");
//   printf("# i.e. both species are products in the same reversible reaction\n");
//   printf("#\n#                  numer   denom\n");
//   printf("# J term   rxn ID  spc ID  spc ID  sparse mtx ID\n");
//   for(j=0; j<(w->termList->nRevDestroy); j++)
//     {
//       printf("%6d  %6d  %6d  %6d  %6d\n",j,
// 	     w->termList->revDestroy[j].rxnIdx,
// 	     w->mtxRowIdx[w->termList->revDestroy[j].sparseIdx],
// 	     w->termList->revDestroy[j].concIdx,
// 	     w->termList->revDestroy[j].sparseIdx);
//     }
//   printf("\n\n");
//   printf("# Reverse Creation Terms:\n");
//   printf("# i.e. numerator species is a reactant and the denominator\n");
//   printf("#      is a product in the same reversible reaction\n");
//   printf("#\n#                  numer   denom\n");
//   printf("# J term   rxn ID  spc ID  spc ID\n");
//   for(j=0; j<(w->termList->nRevCreate); j++)
//     {
//       printf("%6d  %6d  %6d  %6d  %6d\n",j,
// 	     w->termList->revCreate[j].rxnIdx,
// 	     w->mtxRowIdx[w->termList->revCreate[j].sparseIdx],
// 	     w->termList->revCreate[j].concIdx,
// 	     w->termList->revCreate[j].sparseIdx);
//     }
// }

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



void transposeSparseMtxIndices(Jsparse *w)
{
   int j,k,m;
   int nSize = w->nSize;
   int nReactors = w->nReactors;

   //Now Transpose matrix pointers
   // This means we don't have to transpose RHS on backsolve
   int * blockedColSum = new int[w->nSize*nReactors + 1];
   int * blockedRowIdx = new int[w->nNonZero*nReactors];
   blockedColSum[0] = 0;
   blockedColSum[1] = 0;
   for(j = 0; j < w->nSize; ++j)
   {
     for(k = 0; k < nReactors; ++k)
     {
       int nElems = w->mtxColSum[j+1] - w->mtxColSum[j];
       blockedColSum[j*nReactors + k + 1] = blockedColSum[j*nReactors+k] + nElems;
       int origColStart = w->mtxColSum[j];
       int newColStart = blockedColSum[j*nReactors+k];
       for(m = 0; m < nElems; ++m)
       {
          blockedRowIdx[newColStart+m] = w->mtxRowIdx[origColStart+m]*nReactors+k;
       }
     }
   }
   memcpy(w->mtxColSum,blockedColSum,sizeof(int)*(nReactors*w->nSize + 1));
   memcpy(w->mtxRowIdx,blockedRowIdx,sizeof(int)*(nReactors*w->nNonZero));

   memcpy(w->reduceColSum,blockedColSum,sizeof(int)*(nReactors*w->nSize + 1));
   memcpy(w->reduceRowIdx,blockedRowIdx,sizeof(int)*(nReactors*w->nNonZero));

   delete[] blockedColSum;
   delete[] blockedRowIdx;

   //Update special sparse pointers for reactons and temperature row
   // to make Jacobian calculation simpler.
   for(j=0; j<w->termList->nFwdDestroy; ++j)
   {
      int sparseId = w->termList->fwdDestroy[j].sparseIdx;

      int colId = w->termList->fwdDestroy[j].concIdx;
      int currColSum = w->mtxColSum[colId*nReactors+1] - w->mtxColSum[colId*nReactors];
      int colElem = sparseId - w->mtxColSum[colId*nReactors]/nReactors;
      w->termList->fwdDestroy[j].sparseIdx = colElem + w->mtxColSum[colId*nReactors];
   }

   for(j=0; j<w->termList->nFwdCreate; j++)
   {
      int sparseId = w->termList->fwdCreate[j].sparseIdx;

      int colId = w->termList->fwdCreate[j].concIdx;
      int currColSum = w->mtxColSum[colId*nReactors+1] - w->mtxColSum[colId*nReactors];
      int colElem = sparseId - w->mtxColSum[colId*nReactors]/nReactors;
      w->termList->fwdCreate[j].sparseIdx = colElem + w->mtxColSum[colId*nReactors];
   }

   for(j=0; j<w->nSize; j++) // column number
   {
       int sparseId = w->lastRowIdx[j];

       int colElem = sparseId - w->mtxColSum[j*nReactors]/nReactors;
       w->lastRowIdx[j] = colElem + w->mtxColSum[j*nReactors];
   }

   for(j=0; j<w->nSize; j++) // column number
   {
       int sparseId = w->diagIdx[j];

       int colElem = sparseId - w->mtxColSum[j*nReactors]/nReactors;
       w->diagIdx[j] = colElem + w->mtxColSum[j*nReactors];
   }
}
