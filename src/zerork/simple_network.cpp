#include <stdlib.h>
#include "simple_network.h"
#include "constants.h"

namespace zerork {

#ifdef EXIT_THROWS_EXCEPTION
  // create a local function to overide the system exit and throw an exception
  // with the status integer.
  static void exit(int status) {throw status;}
#endif // EXIT_THROWS_EXCEPTION

simple_network::simple_network(ckr::CKReader *ckrobj,
				       nasa_poly_group &tPtr)
{
  if(ckrobj == NULL)
    {
      printf("ERROR: in simple_network() CKReader object not allocated\n");
      exit(-1);
    }
  if(&tPtr == NULL)
    {
      printf("ERROR: in simple_network() nasa_poly_group object not allocated\n");
      exit(-1);
    }

  // determine the multiplier to convert the activation energy units in the
  // prescribed mech into kelvin
  if(ckrobj->units.ActEnergy == ckr::Cal_per_Mole)
    {convertE = CAL_PER_MOL_TACT;}
  else if(ckrobj->units.ActEnergy == ckr::Kelvin)
    {convertE = 1.0;}
  else if(ckrobj->units.ActEnergy == ckr::Kjoules_per_Mole)
    {convertE = KJOULES_PER_MOL_TACT;}
  else
    {
      printf("ERROR: mechanism Activation Energy units type %d not recognized\n",
	     ckrobj->units.ActEnergy);
      exit(-1);
    }

  // determine the multiplier to convert the concentration units in the rate
  // constants to kmol/m^3
  if(ckrobj->units.Quantity == ckr::Moles)
    {
      // note that the length scale in the CKReader object is [cm]
      // [mol/cm^3] * convertC = [kmol/m^3]
      convertC=1000.0;
    }
   else
    {
      printf("ERROR: mechanism Quantity units type %d not recognized\n",
	     ckrobj->units.Quantity);
      exit(-1);
    }

  thermoPtr=&tPtr;

  setRxnStepIdxMaps(ckrobj);
  setParticipants(ckrobj);

}

simple_network::~simple_network()
{
  delete [] productSpcIdx;
  delete [] reactantSpcIdx;
  delete [] addrProductOfStep;
  delete [] addrReactantOfStep;
  delete [] nProductOfStep;
  delete [] nReactantOfStep;
  delete [] rxnIdxOfStep;
  delete [] stepIdxOfRevRxn;
  delete [] stepIdxOfFwdRxn;
}


int simple_network::setRxnStepIdxMaps(ckr::CKReader *ckrobj)
{
  int j,k;

  nRxn=static_cast<int>(ckrobj->reactions.size());

  stepIdxOfFwdRxn = new int[nRxn];
  if(stepIdxOfFwdRxn == NULL)
    {
      printf("ERROR: allocation failure in simple_network::setRxnStepIdxMaps()\n");
      return -1;
    }
  stepIdxOfRevRxn = new int[nRxn];
  if(stepIdxOfRevRxn == NULL)
    {
      printf("ERROR: allocation failure in simple_network::setRxnStepIdxMaps()\n");
      return -2;
    }

  k=0;
  for(j=0; j<nRxn; j++)
    {
      stepIdxOfFwdRxn[j] = k; ++k;

      stepIdxOfRevRxn[j]=MIN_INT32;
      if(ckrobj->reactions[j].isReversible &&
	 ckrobj->reactions[j].krev.A != 0.0)
	{stepIdxOfRevRxn[j]=k; ++k;}
    }
  nStep=k;

  rxnIdxOfStep = new int[nStep];
  if(rxnIdxOfStep == NULL)
    {
      printf("ERROR: allocation failure in simple_network::setRxnStepIdxMaps()\n");
      return -3;
    }

  k=0;
  for(j=0; j<nRxn; j++)
    {
      rxnIdxOfStep[k] = j; ++k;

      if(ckrobj->reactions[j].isReversible &&
	 ckrobj->reactions[j].krev.A != 0.0)
	{rxnIdxOfStep[k]=j+nRxn; ++k;}
    }
  return 0;
}

int simple_network::setParticipants(ckr::CKReader *ckrobj)
{
  int j,k,m,rxnIdx,spcIdx,spcCtr,intCount;
  double reactantCount, productCount, count;
  string spcName;


  nReactantOfStep    = new int[nStep];
  if(nReactantOfStep == NULL)
    {
      printf("ERROR: allocation failure in simple_network::setParticipants()\n");
      return -1;
    }

  nProductOfStep     = new int[nStep];
  if(nProductOfStep == NULL)
    {
      printf("ERROR: allocation failure in simple_network::setParticipants()\n");
      return -2;
    }
  addrReactantOfStep = new int[nStep];
  if(addrReactantOfStep == NULL)
    {
      printf("ERROR: allocation failure in simple_network::setParticipants()\n");
      return -3;
    }

  addrProductOfStep  = new int[nStep];
  if(addrProductOfStep  == NULL)
    {
      printf("ERROR: allocation failure in simple_network::setParticipants()\n");
      return -4;
    }

  totProduct=totReactant=0;
  for(j=0; j<nStep; j++)
    {
      rxnIdx=rxnIdxOfStep[j];

      if(rxnIdx < nRxn) // forward reaction direction
	{
	  // scan the reactant and product counts
	  reactantCount=0.0;
	  for(k=0; k<(int)ckrobj->reactions[rxnIdx].reactants.size(); k++)
	    {
	      count=ckrobj->reactions[rxnIdx].reactants[k].number;
	      intCount=(int)(count+0.5);
	      if(fabs(count-(double)intCount) >= STOICH_TOL)
		{
		  printf("ERROR: non-integer stoichiometric coefficient\n");
		  printf("       found in reactant %d of fwd reaction %d\n",
			 k,rxnIdx);
		  exit(-1);
		}
	      reactantCount+=count;
	    }

	  productCount=0.0;
	  for(k=0; k<(int)ckrobj->reactions[rxnIdx].products.size(); k++)
	    {
	      count=ckrobj->reactions[rxnIdx].products[k].number;
	      intCount=(int)(count+0.5);
	      if(fabs(count-(double)intCount) >= STOICH_TOL)
		{
		  printf("ERROR: non-integer stoichiometric coefficient\n");
		  printf("       found in product %d of fwd reaction %d\n",
			 k,rxnIdx);
		  exit(-1);
		}
	      productCount+=count;
	    }
	}
      else // reverse reaction direction
	{
	  rxnIdx-=nRxn;
	  // scan the reactant and product counts
	  reactantCount=0.0;
	  for(k=0; k<(int)ckrobj->reactions[rxnIdx].products.size(); k++)
	    {reactantCount+=ckrobj->reactions[rxnIdx].products[k].number;}

	  productCount=0.0;
	  for(k=0; k<(int)ckrobj->reactions[rxnIdx].reactants.size(); k++)
	    {productCount+=ckrobj->reactions[rxnIdx].reactants[k].number;}
	}

      nReactantOfStep[j] = (int)(reactantCount + 0.5);
      totReactant+=nReactantOfStep[j];

      nProductOfStep[j]  = (int)(productCount  + 0.5);
      totProduct+=nProductOfStep[j];
    }

  addrReactantOfStep[0]=0;
  addrProductOfStep[0]=0;
  maxReactantInStep=nReactantOfStep[0];
  maxProductInStep=nProductOfStep[0];

  for(j=1; j<nStep; j++)
    {
      addrReactantOfStep[j]=addrReactantOfStep[j-1]+nReactantOfStep[j-1];
      if(nReactantOfStep[j] > maxReactantInStep)
	{maxReactantInStep = nReactantOfStep[j];}

      addrProductOfStep[j]=addrProductOfStep[j-1]+nProductOfStep[j-1];
      if(nProductOfStep[j] > maxProductInStep)
	{maxProductInStep = nProductOfStep[j];}
    }

  reactantSpcIdx = new int[totReactant];

  if(reactantSpcIdx == NULL)
    {
      printf("ERROR: allocation failure in simple_network::setParticipants()\n");
      return -5;
    }
  productSpcIdx = new int[totProduct];
  if(productSpcIdx == NULL)
    {
      printf("ERROR: allocation failure in simple_network::setParticipants()\n");
      return -6;
    }

  for(j=0; j<nStep; j++)
    {

      rxnIdx=rxnIdxOfStep[j];

      if(rxnIdx < nRxn) // forward reaction direction
	{
	  spcCtr=0;
	  for(k=0; k<(int)ckrobj->reactions[rxnIdx].reactants.size(); k++)
	    {
	      spcName=ckrobj->reactions[rxnIdx].reactants[k].name;
	      spcIdx=ckrobj->speciesData[spcName].index;
	      intCount=
		(int)(ckrobj->reactions[rxnIdx].reactants[k].number+0.5);
	      for(m=0; m<intCount; m++)
		{
		  reactantSpcIdx[addrReactantOfStep[j]+spcCtr]=spcIdx;
		  ++spcCtr;
		}
	    }
	  if(spcCtr != nReactantOfStep[j])
	    {
	      printf("ERROR: the reactant species count in step %d does not\n",
		     j);
	      printf("       match the initial count (init %d != %d)\n",
		     nReactantOfStep[j],spcCtr);
	      return -7;
	    }

	  spcCtr=0;
	  for(k=0; k<(int)ckrobj->reactions[rxnIdx].products.size(); k++)
	    {
	      spcName=ckrobj->reactions[rxnIdx].products[k].name;
	      spcIdx=ckrobj->speciesData[spcName].index;
	      intCount=
		(int)(ckrobj->reactions[rxnIdx].products[k].number+0.5);
	      for(m=0; m<intCount; m++)
		{
		  productSpcIdx[addrProductOfStep[j]+spcCtr]=spcIdx;
		  ++spcCtr;
		}
	    }
	  if(spcCtr != nProductOfStep[j])
	    {
	      printf("ERROR: the product species count in step %d does not\n",
		     j);
	      printf("       match the initial count (init %d != %d)\n",
		     nProductOfStep[j],spcCtr);
	      return -8;
	    }
 	}
      else // reverse reaction direction
	{
	  rxnIdx-=nRxn;

	  spcCtr=0;
	  for(k=0; k<(int)ckrobj->reactions[rxnIdx].products.size(); k++)
	    {
	      spcName=ckrobj->reactions[rxnIdx].products[k].name;
	      spcIdx=ckrobj->speciesData[spcName].index;
	      intCount=
		(int)(ckrobj->reactions[rxnIdx].products[k].number+0.5);
	      for(m=0; m<intCount; m++)
		{
		  reactantSpcIdx[addrReactantOfStep[j]+spcCtr]=spcIdx;
		  ++spcCtr;
		}
	    }
	  spcCtr=0;
	  for(k=0; k<(int)ckrobj->reactions[rxnIdx].reactants.size(); k++)
	    {
	      spcName=ckrobj->reactions[rxnIdx].reactants[k].name;
	      spcIdx=ckrobj->speciesData[spcName].index;
	      intCount=
		(int)(ckrobj->reactions[rxnIdx].reactants[k].number+0.5);
	      for(m=0; m<intCount; m++)
		{
		  productSpcIdx[addrProductOfStep[j]+spcCtr]=spcIdx;
		  ++spcCtr;
		}
	    }
	}
   }

  return 0;
}

int setTemperatureTypeMaps(ckr::CKReader *ckrobj);
{
  int *fwdArrheniusDistinct;
  int *revArrheniusDistinct;
  int *auxArrheniusDistinct;
  double *A;
  double *Texp;
  double *Tact;

  fwdArrheniusDistinct=new int[nRxn];
  revArrheniusDistinct=new int[nRxn];
  auxArrheniusDistinct=new int[nRxn];
  A    = new double[nRxn*3];
  Texp = new double[nRxn*3];
  Tact =
}



void simple_network::print_info(const species *spcList)
{
  int j,k,spcIdx;
  printf("# From void simple_network::print_info(const species *)\n");
  printf("#\n");
  printf("# activation energy conversion const to    [K] : %.18g\n",convertE);
  printf("# concentration conversion const to [kmol/m^3] : %.18g\n",convertC);
  printf("# number of chemical reactions                 : %d\n",nRxn);
  printf("# number of chemical reaction steps (1-way)    : %d\n",nStep);
  printf("# total number of reactants                    : %d\n",totReactant);
  printf("# total number of products                     : %d\n",totProduct);
  printf("# max   number of reactants in a step          : %d\n",
	 maxReactantInStep);
  printf("# max   number of products  in a step          : %d\n",
	 maxProductInStep);
  printf("# Reaction steps:\n");

  if(spcList == NULL)
    {
      for(j=0; j<nStep; j++)
	{
	  printf("#   step[%d] (rxn %d ",j,rxnIdxOfStep[j]%nRxn);
	  if(rxnIdxOfStep[j]<nRxn)
	    {printf("fwd): ");}
	  else
	    {printf("rev): ");}

	  spcIdx=reactantSpcIdx[addrReactantOfStep[j]];

	  printf("%d ",spcIdx);
	  for(k=1; k<nReactantOfStep[j]; k++)
	    {
	      spcIdx=reactantSpcIdx[addrReactantOfStep[j]+k];
	      printf("+ %d ",spcIdx);}

	  spcIdx=productSpcIdx[addrProductOfStep[j]];
	  printf("-> %d ",spcIdx);
	  for(k=1; k<nProductOfStep[j]; k++)
	    {
	      spcIdx=productSpcIdx[addrProductOfStep[j]+k];
	      printf("+ %d ",spcIdx);
	    }
	  printf("\n");
	}
    }
  else
    {
     for(j=0; j<nStep; j++)
	{
	  printf("#   step[%d] (rxn %d ",j,rxnIdxOfStep[j]%nRxn);
	  if(rxnIdxOfStep[j]<nRxn)
	    {printf("fwd): ");}
	  else
	    {printf("rev): ");}

	  spcIdx=reactantSpcIdx[addrReactantOfStep[j]];

	  printf("%s ",spcList[spcIdx].getName_c_str());
	  for(k=1; k<nReactantOfStep[j]; k++)
	    {
	      spcIdx=reactantSpcIdx[addrReactantOfStep[j]+k];
	      printf("+ %s ",spcList[spcIdx].getName_c_str());}

	  spcIdx=productSpcIdx[addrProductOfStep[j]];
	  printf("-> %s ",spcList[spcIdx].getName_c_str());
	  for(k=1; k<nProductOfStep[j]; k++)
	    {
	      spcIdx=productSpcIdx[addrProductOfStep[j]+k];
	      printf("+ %s ",spcList[spcIdx].getName_c_str());
	    }
	  printf("\n");
	}
    }


}

} // namespace zerork
