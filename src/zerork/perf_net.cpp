#include <stdlib.h>
#include <assert.h>

#include "utilities.h"
#include "perf_net.h"


namespace zerork {

static void UnsupportedFeature(const char filename[], const int line_num)
{
  printf("# ERROR: Zero-RK reached an unsupported feature\n");
  printf("#            in file = %s\n",filename);
  printf("#            at line = %d\n",line_num);
  fflush(stdout);
  assert(0);  // disabled with macro #define NDEBUG
              // or -DNDEBUG for the compiler
}

perf_net::perf_net(info_net &netobj, rate_const &Kobj)
{
  int j,k,rCtr,pCtr;

  nStep   = netobj.getNumSteps();
  maxReactants = netobj.getMaxReactantsInStep();
  totReac = netobj.getTotalReactants();
  totProd = netobj.getTotalProducts();
  nSpc    = Kobj.getNumSpecies();

  // defines the products and reactants for steps that have all integer
  // stoichiometric coefficients
  reactantSpcIdxList  = new int[totReac];
  reactantStepIdxList = new int[totReac];
  productSpcIdxList   = new int[totProd];
  productStepIdxList  = new int[totProd];

  stepRate       = new double[nStep];
  netSpcProdRate = new double[nSpc];
  spcDestroyRate = new double[nSpc];
  spcCreateRate  = new double[nSpc];

  rateConstPtr=&Kobj;
  infoPtr=&netobj;

  rCtr = 0; // reactant species counter
  pCtr = 0; // product  species counter
  for(j=0; j<nStep; j++)
  {
    // getOrderOfStep(j) returns -1 for non-integer reactions
    for(k=0; k<netobj.getOrderOfStep(j); k++)
    {
      reactantSpcIdxList[rCtr]  = netobj.getSpecIdxOfStepReactant(j,k);
      reactantStepIdxList[rCtr] = j;
      ++rCtr;
    }
     // getNumProductsOfStep(j) returns -1 for non-integer reactions
    for(k=0; k<netobj.getNumProductsOfStep(j); k++)
    {
      productSpcIdxList[pCtr]  = netobj.getSpecIdxOfStepProduct(j,k);
      productStepIdxList[pCtr] = j;
      ++pCtr;
    }
  }
  assert(pCtr == totProd);
  assert(rCtr == totReac);


  non_integer_network_ = netobj.getNonIntegerReactionNetwork();
  niReactantSpcIdxList.clear();
  niReactantStepIdxList.clear();
  niReactantStoichNumList.clear();
  niProductSpcIdxList.clear();
  niProductStepIdxList.clear();
  niProductStoichNumList.clear();
  use_non_integer_network_ = true;
  if(non_integer_network_.GetNumNonIntegerSteps() == 0) {
    use_non_integer_network_ = false;
  } else {
    for(j=0; j<nStep; j++)
    {
      if(non_integer_network_.HasStep(j))
      {
        std::vector<int> step_reactant_indexes = non_integer_network_.GetReactantIndexesOfStep(j);
        std::vector<double> step_reactant_stoich_nums = non_integer_network_.GetReactantStoichNumsOfStep(j);
        const int num_reactants = step_reactant_indexes.size();
        if(num_reactants > maxReactants) maxReactants = num_reactants;
        for(k=0; k<num_reactants; k++) 
        {
          niReactantSpcIdxList.push_back(step_reactant_indexes[k]);
          niReactantStoichNumList.push_back(step_reactant_stoich_nums[k]);
          niReactantStepIdxList.push_back(j);
        }
        std::vector<int> step_product_indexes = non_integer_network_.GetProductIndexesOfStep(j);
        std::vector<double> step_product_stoich_nums = non_integer_network_.GetProductStoichNumsOfStep(j);
        const int num_products = step_product_indexes.size();
        for(k=0; k<num_products; k++)
        {
          niProductSpcIdxList.push_back(step_product_indexes[k]);
          niProductStoichNumList.push_back(step_product_stoich_nums[k]);
          niProductStepIdxList.push_back(j);
        }
      }
    }
  }
  niTotProd = niProductSpcIdxList.size();
  niTotReac = niReactantSpcIdxList.size();;

  //Initialize timing data
  cpuKTime = cpuStepTime = cpuProdTime = cpuNetTime = 0.0;

  use_external_rates = false;
}

perf_net::~perf_net()
{
  delete [] spcCreateRate;
  delete [] spcDestroyRate;
  delete [] netSpcProdRate;
  delete [] stepRate;
  delete [] productStepIdxList;
  delete [] productSpcIdxList;
  delete [] reactantStepIdxList;
  delete [] reactantSpcIdxList;
}

void perf_net::calcRatesFromTC(const double T, const double C[],
                                   double netOut[], double createOut[],
                                   double destroyOut[], double stepOut[])
{
  int j,loopLim; //,reacId,stepId;
  double startTime,tmpTime;

  startTime = getHighResolutionTime();
  rateConstPtr->updateK(T,&C[0],&stepOut[0]); // store K(T,p,C) in stepOut[]
  cpuKTime += getHighResolutionTime() - startTime;

  if(use_external_rates)
  {
      UnsupportedFeature(__FILE__, __LINE__);

      startTime = getHighResolutionTime();
      (*ex_func_calc_rates)(&C[0],&stepOut[0],&createOut[0],&destroyOut[0]);
      tmpTime = getHighResolutionTime() - startTime;
      //FIXME: Poor approximation below.
      cpuStepTime += tmpTime/2;
      cpuProdTime += tmpTime/2;
  }
  else
  {
      startTime = getHighResolutionTime();
      // compute the rate of progress of each step
      for(j=0; j<totReac; ++j)
        {stepOut[reactantStepIdxList[j]]*=C[reactantSpcIdxList[j]];}
      cpuStepTime += getHighResolutionTime() - startTime;

      startTime = getHighResolutionTime();
    // set the species creation and destruction rates to zero
//      for(j=0; j<nSpc; ++j)
//        {createOut[j]=destroyOut[j]=0.0;}

      memset(createOut,0,nSpc*sizeof(double));
      memset(destroyOut,0,nSpc*sizeof(double));

      if(use_non_integer_network_) {
        non_integer_network_.UpdateRatesOfProgress(C,stepOut);
	non_integer_network_.GetCreationRates(stepOut,createOut);
        non_integer_network_.GetDestructionRates(stepOut,destroyOut);
      }

      // compute the species destruction rate by adding each steps rate of progress
      // to the sum for each reactant species found
      for(j=0; j<totReac; ++j)
        {destroyOut[reactantSpcIdxList[j]]+=stepOut[reactantStepIdxList[j]];}

      // compute the species creation rate by adding each steps rate of progress
      // to the sum for each product species found
      for(j=0; j<totProd; ++j)
        {createOut[productSpcIdxList[j]]+=stepOut[productStepIdxList[j]];}
      cpuProdTime += getHighResolutionTime() - startTime;
  }

  startTime = getHighResolutionTime();
  // compute the net species production rate = create - destroy
  loopLim = nSpc;
  for(j=0; j<loopLim; ++j)
    {netOut[j]=createOut[j]-destroyOut[j];}
  cpuNetTime += getHighResolutionTime() - startTime;
}

void perf_net::calcRatesFromTC_StepLimiter(const double T,
                                           const double C[],
		                           const double step_limiter[],
			                   double netOut[],
                                           double createOut[],
			                   double destroyOut[],
                                           double stepOut[])

{
  int j,loopLim; //,reacId,stepId;
  double startTime,tmpTime;

  startTime = getHighResolutionTime();
  rateConstPtr->updateK(T,&C[0],&stepOut[0]); // store K(T,p,C) in stepOut[]
  cpuKTime += getHighResolutionTime() - startTime;

  // Apply the step_limiter to each K
  const int const_num_steps = nStep;
  for(j=0; j<const_num_steps; ++j) {
    if(0.0 < step_limiter[j] && step_limiter[j] < 1.0e+300) {

      //if(stepOut[j] > step_limiter[j]) {
      //  printf("# INFO: In calcRatesFromTC_StepLimiter(...),\n");
      //  printf("# INFO:   step [%4d] %14.7e exceeds limit of %14.7e\n",
      //         j,stepOut[j],step_limiter[j]);
      //}
      stepOut[j] *= step_limiter[j]/(step_limiter[j]+stepOut[j]);

    }
  }
  fflush(stdout);


  if(use_external_rates)
  {
      UnsupportedFeature(__FILE__, __LINE__);

      startTime = getHighResolutionTime();
      (*ex_func_calc_rates)(&C[0],&stepOut[0],&createOut[0],&destroyOut[0]);
      tmpTime = getHighResolutionTime() - startTime;
      //FIXME: Poor approximation below.
      cpuStepTime += tmpTime/2;
      cpuProdTime += tmpTime/2;
  }
  else
  {
      startTime = getHighResolutionTime();
      // compute the rate of progress of each step
      for(j=0; j<totReac; ++j)
        {stepOut[reactantStepIdxList[j]]*=C[reactantSpcIdxList[j]];}
      cpuStepTime += getHighResolutionTime() - startTime;

      startTime = getHighResolutionTime();
    // set the species creation and destruction rates to zero
//      for(j=0; j<nSpc; ++j)
//        {createOut[j]=destroyOut[j]=0.0;}

      memset(createOut,0,nSpc*sizeof(double));
      memset(destroyOut,0,nSpc*sizeof(double));

      if(use_non_integer_network_) {
        non_integer_network_.UpdateRatesOfProgress(C,stepOut);
	non_integer_network_.GetCreationRates(stepOut,createOut);
        non_integer_network_.GetDestructionRates(stepOut,destroyOut);
      }

      // compute the species destruction rate by adding each steps rate of progress
      // to the sum for each reactant species found
      for(j=0; j<totReac; ++j)
        {destroyOut[reactantSpcIdxList[j]]+=stepOut[reactantStepIdxList[j]];}

      // compute the species creation rate by adding each steps rate of progress
      // to the sum for each product species found
      for(j=0; j<totProd; ++j)
        {createOut[productSpcIdxList[j]]+=stepOut[productStepIdxList[j]];}
      cpuProdTime += getHighResolutionTime() - startTime;
  }

  startTime = getHighResolutionTime();
  // compute the net species production rate = create - destroy
  loopLim = nSpc;
  for(j=0; j<loopLim; ++j)
    {netOut[j]=createOut[j]-destroyOut[j];}
  cpuNetTime += getHighResolutionTime() - startTime;
}

void perf_net::calcRatesFromTC_StepLimiter_perturbROP(const double T,
                                                      const double C[],
                                                      const double step_limiter[],
                                                      const double perturbMult[],
                                                      double netOut[],
                                                      double createOut[],
                                                      double destroyOut[],
                                                      double stepOut[])

{
  int j,loopLim; //,reacId,stepId;
  double startTime,tmpTime;

  startTime = getHighResolutionTime();
  rateConstPtr->updateK(T,&C[0],&stepOut[0]); // store K(T,p,C) in stepOut[]
  cpuKTime += getHighResolutionTime() - startTime;

  // Apply the step_limiter to each K
  const int const_num_steps = nStep;
  for(j=0; j<const_num_steps; ++j) {
    if(0.0 < step_limiter[j] && step_limiter[j] < 1.0e+300) {

      //if(stepOut[j] > step_limiter[j]) {
      //  printf("# INFO: In calcRatesFromTC_StepLimiter(...),\n");
      //  printf("# INFO:   step [%4d] %14.7e exceeds limit of %14.7e\n",
      //         j,stepOut[j],step_limiter[j]);
      //}
      stepOut[j] *= step_limiter[j]/(step_limiter[j]+stepOut[j]);

    }
  }
  fflush(stdout);

  // apply the multiplicative perturbation to the ROP array
  for(j=0; j<nStep; j++)
    {stepOut[j]*=perturbMult[j];}


  if(use_external_rates)
  {
      UnsupportedFeature(__FILE__, __LINE__);

      startTime = getHighResolutionTime();
      (*ex_func_calc_rates)(&C[0],&stepOut[0],&createOut[0],&destroyOut[0]);
      tmpTime = getHighResolutionTime() - startTime;
      //FIXME: Poor approximation below.
      cpuStepTime += tmpTime/2;
      cpuProdTime += tmpTime/2;
  }
  else
  {
      startTime = getHighResolutionTime();
      // compute the rate of progress of each step
      for(j=0; j<totReac; ++j)
        {stepOut[reactantStepIdxList[j]]*=C[reactantSpcIdxList[j]];}
      cpuStepTime += getHighResolutionTime() - startTime;

      startTime = getHighResolutionTime();
    // set the species creation and destruction rates to zero
//      for(j=0; j<nSpc; ++j)
//        {createOut[j]=destroyOut[j]=0.0;}

      memset(createOut,0,nSpc*sizeof(double));
      memset(destroyOut,0,nSpc*sizeof(double));

      if(use_non_integer_network_) {
        non_integer_network_.UpdateRatesOfProgress(C,stepOut);
	non_integer_network_.GetCreationRates(stepOut,createOut);
        non_integer_network_.GetDestructionRates(stepOut,destroyOut);
      }

      // compute the species destruction rate by adding each steps rate of progress
      // to the sum for each reactant species found
      for(j=0; j<totReac; ++j)
        {destroyOut[reactantSpcIdxList[j]]+=stepOut[reactantStepIdxList[j]];}

      // compute the species creation rate by adding each steps rate of progress
      // to the sum for each product species found
      for(j=0; j<totProd; ++j)
        {createOut[productSpcIdxList[j]]+=stepOut[productStepIdxList[j]];}
      cpuProdTime += getHighResolutionTime() - startTime;
  }

  startTime = getHighResolutionTime();
  // compute the net species production rate = create - destroy
  loopLim = nSpc;
  for(j=0; j<loopLim; ++j)
    {netOut[j]=createOut[j]-destroyOut[j];}
  cpuNetTime += getHighResolutionTime() - startTime;
}

void perf_net::calcRatesFromTCM(const double T,
                                const double C[],
                                const double C_mix,
                                double netOut[],
		                double createOut[],
                                double destroyOut[],
		                double stepOut[])
{
  int j,loopLim; //,reacId,stepId;
  double startTime,tmpTime;

  startTime = getHighResolutionTime();
  rateConstPtr->updateK_TCM(T,
                            &C[0],
                            C_mix,
                            &stepOut[0]); // store K(T,p,C) in stepOut[]
  cpuKTime += getHighResolutionTime() - startTime;

  if(use_external_rates)
  {
      UnsupportedFeature(__FILE__, __LINE__);

      startTime = getHighResolutionTime();
      (*ex_func_calc_rates)(&C[0],&stepOut[0],&createOut[0],&destroyOut[0]);
      tmpTime = getHighResolutionTime() - startTime;
      //FIXME: Poor approximation below.
      cpuStepTime += tmpTime/2;
      cpuProdTime += tmpTime/2;
  }
  else
  {
      startTime = getHighResolutionTime();
      // compute the rate of progress of each step
      for(j=0; j<totReac; ++j)
        {stepOut[reactantStepIdxList[j]]*=C[reactantSpcIdxList[j]];}
      cpuStepTime += getHighResolutionTime() - startTime;

      startTime = getHighResolutionTime();
    // set the species creation and destruction rates to zero
//      for(j=0; j<nSpc; ++j)
//        {createOut[j]=destroyOut[j]=0.0;}
      memset(createOut,0,nSpc*sizeof(double));
      memset(destroyOut,0,nSpc*sizeof(double));

      if(use_non_integer_network_) {
        non_integer_network_.UpdateRatesOfProgress(C,stepOut);
	non_integer_network_.GetCreationRates(stepOut,createOut);
        non_integer_network_.GetDestructionRates(stepOut,destroyOut);
      }

      // compute the species destruction rate by adding each steps rate of progress
      // to the sum for each reactant species found
      for(j=0; j<totReac; ++j)
        {destroyOut[reactantSpcIdxList[j]]+=stepOut[reactantStepIdxList[j]];}

      // compute the species creation rate by adding each steps rate of progress
      // to the sum for each product species found
      for(j=0; j<totProd; ++j)
        {createOut[productSpcIdxList[j]]+=stepOut[productStepIdxList[j]];}
      cpuProdTime += getHighResolutionTime() - startTime;
  }

  startTime = getHighResolutionTime();
  // compute the net species production rate = create - destroy
  loopLim = nSpc;
  for(j=0; j<loopLim; ++j)
    {netOut[j]=createOut[j]-destroyOut[j];}
  cpuNetTime += getHighResolutionTime() - startTime;
}


void perf_net::calcRatesFromExplicit(const double T, const double C[],
                                   double netOut[], double createOut[],
                                   double destroyOut[], double stepOut[])
{
  UnsupportedFeature(__FILE__, __LINE__);
  int j; //,reacId,stepId;
  rateConstPtr->updateK(T,&C[0],&stepOut[0]); // store K(T,p,C) in stepOut[]
//  rateConstPtr->updateKExplicit(T,&C[0],&stepOut[0]); // store K(T,p,C) in stepOut[]

  (*ex_func_calc_rates)(&C[0],&stepOut[0],&createOut[0],&destroyOut[0]);

  // compute the net species production rate = create - destroy
  for(j=0; j<nSpc; ++j)
    {netOut[j]=createOut[j]-destroyOut[j];}

}

void perf_net::calcRatesFromTC_perturbROP(const double T, const double C[],
           const double perturbMult[], double netOut[], double createOut[],
           double destroyOut[], double stepOut[])
{
  int j; //,reacId,stepId;
  rateConstPtr->updateK(T,&C[0],&stepOut[0]); // store K(T,p,C) in stepOut[]

  // compute the rate of progress of each step
  for(j=0; j<totReac; ++j)
    {stepOut[reactantStepIdxList[j]]*=C[reactantSpcIdxList[j]];}


  // apply the multiplicative perturbation to the ROP array
  for(j=0; j<nStep; j++)
    {stepOut[j]*=perturbMult[j];}

  // set the species creation and destruction rates to zero
  for(j=0; j<nSpc; ++j)
    {createOut[j]=destroyOut[j]=0.0;}

  if(use_non_integer_network_) {
    non_integer_network_.UpdateRatesOfProgress(C,stepOut);
    non_integer_network_.GetCreationRates(stepOut,createOut);
    non_integer_network_.GetDestructionRates(stepOut,destroyOut);
  }

  // compute the species destruction rate by adding each steps rate of progress
  // to the sum for each reactant species found
  for(j=0; j<totReac; ++j)
    {destroyOut[reactantSpcIdxList[j]]+=stepOut[reactantStepIdxList[j]];}

  // compute the species creation rate by adding each steps rate of progress
  // to the sum for each product species found
  for(j=0; j<totProd; ++j)
    {createOut[productSpcIdxList[j]]+=stepOut[productStepIdxList[j]];}

  // compute the net species production rate = create - destroy
  for(j=0; j<nSpc; ++j)
    {netOut[j]=createOut[j]-destroyOut[j];}

}


void perf_net::calcRatesFromTC_unroll16(const double T, const double C[],
				   double netOut[], double createOut[],
				   double destroyOut[], double stepOut[])
{
  UnsupportedFeature(__FILE__, __LINE__);
  int j,k;
  int loopMax,loopRem;
  rateConstPtr->updateK(T,&C[0],&stepOut[0]); // store K(T,p,C) in stepOut[]

  // compute the rate of progress of each step
  loopMax=totReac/16;
  loopRem=totReac%16;
  k=0;
  for(j=0; j<loopRem; j++)
    {
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
    }
  for(j=0; j<loopMax; j++)
    {
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
      stepOut[reactantStepIdxList[k]]*=C[reactantSpcIdxList[k]]; ++k;
    }



  // set the species creation and destruction rates to zero
  for(j=0; j<nSpc; )
    {createOut[j]=destroyOut[j]=0.0; ++j;}

  // compute the species destruction rate by adding each steps rate of progress
  // to the sum for each reactant species found
  for(j=0; j<totReac; )
    {destroyOut[reactantSpcIdxList[j]]+=stepOut[reactantStepIdxList[j]]; ++j;}

  // compute the species creation rate by adding each steps rate of progress
  // to the sum for each product species found
  for(j=0; j<totProd; )
    {createOut[productSpcIdxList[j]]+=stepOut[productStepIdxList[j]]; ++j;}

  // compute the net species production rate = create - destroy
  for(j=0; j<nSpc; )
    {netOut[j]=createOut[j]-destroyOut[j]; ++j;}

}

void perf_net::writeExternalFuncs()
{
  UnsupportedFeature(__FILE__, __LINE__);

  FILE *fptr=fopen("external_funcs.cpp","w");
  if(fptr==NULL) {
      printf("ERROR: could not open file external_funcs.cpp for write operation\n");
      printf("       in perf_net::writeExternalFuncs(...)\n");
      return;
  }

  fprintf(fptr,"#include \"external_funcs.h\"\n");
  fprintf(fptr,"#include \"stdio.h\"\n");
  fprintf(fptr,"#include \"math.h\"\n");
  fprintf(fptr,"#include \"fast_exps.h\"\n");
  fprintf(fptr,"\n\n");

  write_func_check(fptr);
  fprintf(fptr,"\n\n");
  write_func_rates(fptr);
  fprintf(fptr,"\n\n");
  rateConstPtr->write_funcs(fptr);

  fclose(fptr);
}

void perf_net::write_func_check(FILE* fptr)
{
  UnsupportedFeature(__FILE__, __LINE__);

  fprintf(fptr,"bool external_func_check(const int nsp, const int nstep)\n");
  fprintf(fptr,"{\n");

  fprintf(fptr,"  if(nsp != %d || nstep != %d)\n",nSpc,nStep);
  fprintf(fptr,"  {\n");
  fprintf(fptr,"      printf(\"Current mechanism doesn't match external lib. Quitting.\\n\");\n");
  fprintf(fptr,"      return false;\n");
  fprintf(fptr,"  }\n");
  fprintf(fptr,"  return true;\n");
  fprintf(fptr,"}\n");

  fprintf(fptr,"}\n");
}


void perf_net::write_func_rates(FILE* fptr)
{
  UnsupportedFeature(__FILE__, __LINE__);

  int j,k,m,tmpIdx;
  char ropVarName[]="phi";
  char concVarName[]="C";
  char createVarName[]="cdot";
  char destroyVarName[]="ddot";

  int *stepReactantList;
  vector< vector<int> > cdotRoPlist(nSpc, vector <int>(1,-1));
  vector< vector<int> > ddotRoPlist(nSpc, vector <int>(1,-1));

  fprintf(fptr,"void external_func_rates(const double %s[], double %s[], double %s[], double %s[])\n",
                concVarName,ropVarName,createVarName,destroyVarName);
  fprintf(fptr,"{\n");
  fprintf(fptr,"  // compute the rate of progress for each step\n");

  stepReactantList=new int[infoPtr->getMaxReactantsInStep()];

  for(j=0; j<nStep; j++)
    {
      for(k=0; k<infoPtr->getOrderOfStep(j); k++)
	{stepReactantList[k] = infoPtr->getSpecIdxOfStepReactant(j,k);}
      // quick bubble sort
      for(k=0; k<infoPtr->getOrderOfStep(j)-1; k++)
	{
	  for(m=k+1; m<infoPtr->getOrderOfStep(j); m++)
	    {
	      if(stepReactantList[m]<stepReactantList[k])
		{
		  tmpIdx=stepReactantList[m];
		  stepReactantList[m]=stepReactantList[k];
		  stepReactantList[k]=tmpIdx;
		}
	    }
	}
      fprintf(fptr,"  %s[%5d]*=%s[%5d]",ropVarName,j,concVarName,
	      stepReactantList[0]);
      for(k=1; k<infoPtr->getOrderOfStep(j); k++)
	{fprintf(fptr,"*%s[%5d]",concVarName,stepReactantList[k]);}
      fprintf(fptr,";\n");
    }
  //printf("here 1\n"); fflush(stdout);
  for(j=0; j<nStep; j++)
    {
      for(k=0; k<infoPtr->getOrderOfStep(j); k++)
	{
	  //printf("reactant j=%d, k=%d\n",j,k);fflush(stdout);
	  tmpIdx=infoPtr->getSpecIdxOfStepReactant(j,k);
	  ddotRoPlist[tmpIdx].push_back(j);
	}
      for(k=0; k<infoPtr->getNumProductsOfStep(j); k++)
	{
	  //printf("product j=%d, k=%d\n",j,k);fflush(stdout);
	  tmpIdx=infoPtr->getSpecIdxOfStepProduct(j,k);
	  cdotRoPlist[tmpIdx].push_back(j);
	}
    }
  //printf("here 2\n"); fflush(stdout);
  fprintf(fptr,"  // compute the total destruction rates for each species\n");
  for(j=0; j<nSpc; j++)
    {
      if(ddotRoPlist[j].size() == 1)
	{fprintf(fptr,"  %s[%5d]=0.0;\n",destroyVarName,j);}
      else
	{
	  fprintf(fptr,"  %s[%5d]=%s[%5d]",destroyVarName,j,ropVarName,
		  ddotRoPlist[j][1]);
	  for(k=2; k< (int) ddotRoPlist[j].size(); k++)
	    {
	      if(k%5==0)
		{fprintf(fptr,"\n             ");}
	      fprintf(fptr,"+%s[%5d]",ropVarName,ddotRoPlist[j][k]);

	    }
	  fprintf(fptr,";\n");
	}
    }
  //printf("here 3\n"); fflush(stdout);
  fprintf(fptr,"  // compute the total creation rates for each species\n");
  for(j=0; j<nSpc; j++)
    {
      if(cdotRoPlist[j].size() == 1)
	{fprintf(fptr,"  %s[%5d]=0.0;\n",createVarName,j);}
      else
	{
	  fprintf(fptr,"  %s[%5d]=%s[%5d]",createVarName,j,ropVarName,
		  cdotRoPlist[j][1]);
	  for(k=2; k< (int) cdotRoPlist[j].size(); k++)
	    {
	      if(k%5==0)
		{fprintf(fptr,"\n             ");}
	      fprintf(fptr,"+%s[%5d]",ropVarName,cdotRoPlist[j][k]);
	    }
	  fprintf(fptr,";\n");
	}
    }

  fprintf(fptr,"}\n");


  delete [] stepReactantList;

};


} // namespace zerork
