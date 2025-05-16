#include <string.h>
#include <assert.h>
#include "info_net.h"
#include "constants.h"

namespace zerork {

info_net::info_net(ckr::CKReader *ckrobj)
{
  assert(ckrobj != NULL);
 
  setRxnStepIdxMaps(ckrobj);
  setReactionFlags(ckrobj);
  setParticipants(ckrobj);
  SetEnhancedSpeciesList(ckrobj);
}

info_net::~info_net()
{
  delete [] isEnhanced;
  delete [] isFalloff;
  delete [] isThirdBody;
  delete [] isReversible;
  delete [] isNonIntegerStoich;
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

double info_net::getRealOrderOfStep(const int step_id) const
{
  const int reaction_id=getRxnIdxOfStep(step_id);
  if(isNonIntegerStoich[reaction_id] == 0) {
    return static_cast<double>(getOrderOfStep(step_id));
  }
  return non_integer_network.GetOrderOfStep(step_id);

}

int info_net::setRxnStepIdxMaps(ckr::CKReader *ckrobj)
{
  int j,k;

  nRxn=static_cast<int>(ckrobj->reactions.size());

  stepIdxOfFwdRxn = new int[nRxn];
  if(stepIdxOfFwdRxn == NULL)
  {
    printf("ERROR: allocation failure in info_net::setRxnStepIdxMaps()\n");
    return -1;
  }  
  stepIdxOfRevRxn = new int[nRxn];
  if(stepIdxOfRevRxn == NULL)
  {
    printf("ERROR: allocation failure in info_net::setRxnStepIdxMaps()\n");
    return -2;
  }

  k=0;
  for(j=0; j<nRxn; j++)
  {
    stepIdxOfFwdRxn[j] = k; ++k;

    stepIdxOfRevRxn[j]=MIN_INT32;
    if(ckrobj->reactions[j].isReversible && ckrobj->reactions[j].krev.A != 0.0)
    {
      stepIdxOfRevRxn[j]=k;
      ++k;
    }
  }
  nStep=k;
  
  rxnIdxOfStep = new int[nStep];
  if(rxnIdxOfStep == NULL)
  {
    printf("ERROR: allocation failure in info_net::setRxnStepIdxMaps()\n");
    return -3;
  }

  k=0;
  for(j=0; j<nRxn; j++)
  {
    rxnIdxOfStep[k] = j; ++k;

    if(ckrobj->reactions[j].isReversible && 
       ckrobj->reactions[j].krev.A != 0.0)
    {
      rxnIdxOfStep[k]=j+nRxn;
      ++k;
    }
  }
  return 0;
}

int info_net::SetEnhancedSpeciesList(ckr::CKReader *ckrobj)
{
  const int num_reactions_const = static_cast<int>(ckrobj->reactions.size());
  const int num_species_const   = static_cast<int>(ckrobj->species.size());
  int num_enhanced_reactions=0;
  int num_enhanced_species;
  string species_name;
  int species_index;
  EnhancedSpecies new_step_enhanced_species;
  std::map<string, double>::iterator reader_iter;
  std::map<string, int> index_of_species_named;
  std::map<string, int>::iterator name_iter;

  // initialize the array of flags indicating reactions with enhanced
  // third body terms
  isEnhanced = new int[num_reactions_const];
  for(int j=0; j<num_reactions_const; ++j) {
    isEnhanced[j]=0;
  }  

  // create index map for species names
  for(int j=0; j<num_species_const; ++j) {
    index_of_species_named[ckrobj->species[j].name] = j;
  }

  enhanced_species_of_step.clear();

  for(int j=0; j<num_reactions_const; ++j) {
    num_enhanced_species = ckrobj->reactions[j].e3b.size();

    new_step_enhanced_species.index.clear();
    new_step_enhanced_species.alpha.clear();

    if(num_enhanced_species > 0) {

      reader_iter = ckrobj->reactions[j].e3b.begin();
      for(int k=0; k<num_enhanced_species; ++k) {

        species_name = (*reader_iter).first;
        name_iter = index_of_species_named.find(species_name);

        if(name_iter == index_of_species_named.end()) {
          // species name not found
	  printf("WARNING: reaction %d enhanced third body species %s\n",
		 j,species_name.c_str());
	  printf("	   not found in mechanism. Ignoring it.\n");
	} else {
          
          species_index = index_of_species_named[species_name];
          new_step_enhanced_species.index.push_back(species_index);

          // alpha is the 3rd body enchancement factor minus one
          new_step_enhanced_species.alpha.push_back((*reader_iter).second-1.0);
          //printf("# Reaction %d, species %s, enhancement - 1.0: %5.3f\n",
          //       j,species_name.c_str(),(*reader_iter).second-1.0);
        }
        *reader_iter++;     
      }
       isEnhanced[j] = 1; // set enhanced species flag
      ++num_enhanced_reactions;
    } // end if(num_enhanced_species > 0)
    enhanced_species_of_step.push_back(new_step_enhanced_species); 
  }

  return num_enhanced_reactions;
}

int info_net::getEnhancementFactorsOfStep(const int id,
                                          std::vector<int> *species_id,
                                          std::vector<double> *alpha)
{
  species_id->clear();
  alpha->clear();
  *species_id = enhanced_species_of_step[id].index;
  *alpha      = enhanced_species_of_step[id].alpha;
  return enhanced_species_of_step[id].index.size(); 
}

int info_net::setParticipants(ckr::CKReader *ckrobj)
{
  int j,k,m,rxnIdx,spcIdx,spcCtr,intCount;
  double reactantCount, productCount, count;
  string spcName;
 


  nReactantOfStep    = new int[nStep];
  if(nReactantOfStep == NULL)
  {
    printf("ERROR: allocation failure in info_net::setParticipants()\n");
    return -1;
  }

  nProductOfStep     = new int[nStep];
  if(nProductOfStep == NULL)
  {
    printf("ERROR: allocation failure in info_net::setParticipants()\n");
    return -2;
  }
  addrReactantOfStep = new int[nStep];
  if(addrReactantOfStep == NULL)
  {
    printf("ERROR: allocation failure in info_net::setParticipants()\n");
    return -3;
  }

  addrProductOfStep  = new int[nStep];
  if(addrProductOfStep  == NULL)
  {
    printf("ERROR: allocation failure in info_net::setParticipants()\n");
    return -4;
  }
  // create index map for species names
  const int num_species   = static_cast<int>(ckrobj->species.size());
  std::map<string, int> index_of_species_named;
  for(int j=0; j<num_species; ++j) {
    index_of_species_named[ckrobj->species[j].name] = j;
  }

  totProduct=totReactant=0;
  for(j=0; j<nStep; j++)
  {
    rxnIdx=rxnIdxOfStep[j];
    int lookup_id = ((rxnIdx < nRxn) ? rxnIdx : rxnIdx-nRxn);
     
    if(isNonIntegerStoich[lookup_id] == 0) {
 
      if(rxnIdx < nRxn) // forward reaction direction
      {
        // scan the reactant and product counts
        reactantCount=0.0;
        for(k=0; k<(int)ckrobj->reactions[rxnIdx].reactants.size(); k++)
        {
          count=ckrobj->reactions[rxnIdx].reactants[k].number;
          intCount=(int)(count+0.5);
//          if(fabs(count-(double)intCount) >= STOICH_TOL)
//          {
//            printf("# ERROR: non-integer stoichiometric coefficient\n");
//            printf("#        found in reactant %d of fwd reaction %d\n",
//                   k,rxnIdx);
//          }
          assert(fabs(count-(double)intCount) < STOICH_TOL);
          reactantCount+=count;
        }
     
        productCount=0.0;
        for(k=0; k<(int)ckrobj->reactions[rxnIdx].products.size(); k++)
        {
          count=ckrobj->reactions[rxnIdx].products[k].number;
          intCount=(int)(count+0.5);
//          if(fabs(count-(double)intCount) >= STOICH_TOL)
//          {
//            printf("# ERROR: non-integer stoichiometric coefficient\n");
//            printf("#        found in product %d of fwd reaction %d\n",
//                   k,rxnIdx);
//          }
          assert(fabs(count-(double)intCount) < STOICH_TOL);
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
    } else { // end if(isNonIntegerStoich[rxnIdx] == 0)
      // non-integer reactions do not contribute to the reactant and product
      // counts
      nReactantOfStep[j] = 0;
      nProductOfStep[j] = 0;
      ReactionDirection dir = ((lookup_id == rxnIdx) ? FORWARD : REVERSE);
      non_integer_network.AddStep(ckrobj->reactions[lookup_id],
                                  index_of_species_named,
                                  lookup_id, // reaction index
                                  dir,       // reaction direction 
                                  j);        // step index
                                  
    }
  } // end for(j=0; j<nStep; j++)

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
    printf("ERROR: allocation failure in info_net::setParticipants()\n");
    return -5;
  }
  productSpcIdx = new int[totProduct];
  if(productSpcIdx == NULL)
  {
    printf("ERROR: allocation failure in info_net::setParticipants()\n");
    return -6;
  }

  for(j=0; j<nStep; j++)
  {
    rxnIdx=rxnIdxOfStep[j];
    int lookup_id = ((rxnIdx < nRxn) ? rxnIdx : rxnIdx-nRxn);

    if(isNonIntegerStoich[lookup_id] == 0) {

      if(rxnIdx < nRxn) {// forward reaction direction

        spcCtr=0;
        for(k=0; k<(int)ckrobj->reactions[rxnIdx].reactants.size(); k++) {
          spcName=ckrobj->reactions[rxnIdx].reactants[k].name;
          spcIdx=ckrobj->speciesData[spcName].index;
          intCount= (int)(ckrobj->reactions[rxnIdx].reactants[k].number+0.5);
	  for(m=0; m<intCount; m++) {
            reactantSpcIdx[addrReactantOfStep[j]+spcCtr]=spcIdx;
            ++spcCtr;
          }
        }
        if(spcCtr != nReactantOfStep[j]) {
          printf("ERROR: the reactant species count in step %d does not\n",
                 j);
          printf("       match the initial count (init %d != %d)\n",
                 nReactantOfStep[j],spcCtr);
          return -7;
        }

        spcCtr=0;
        for(k=0; k<(int)ckrobj->reactions[rxnIdx].products.size(); k++) {
          spcName=ckrobj->reactions[rxnIdx].products[k].name;
          spcIdx=ckrobj->speciesData[spcName].index;
          intCount= (int)(ckrobj->reactions[rxnIdx].products[k].number+0.5);
	  for(m=0; m<intCount; m++) {
            productSpcIdx[addrProductOfStep[j]+spcCtr]=spcIdx;
            ++spcCtr;
          }
        }
        if(spcCtr != nProductOfStep[j]) {
          printf("ERROR: the product species count in step %d does not\n",
                 j);
          printf("       match the initial count (init %d != %d)\n",
                 nProductOfStep[j],spcCtr);
          return -8;
        }
      }
      else { // reverse reaction direction

        rxnIdx-=nRxn;

        spcCtr=0;
        for(k=0; k<(int)ckrobj->reactions[rxnIdx].products.size(); k++) {
          spcName=ckrobj->reactions[rxnIdx].products[k].name;
          spcIdx=ckrobj->speciesData[spcName].index;
          intCount= (int)(ckrobj->reactions[rxnIdx].products[k].number+0.5);
	  for(m=0; m<intCount; m++) {
            reactantSpcIdx[addrReactantOfStep[j]+spcCtr]=spcIdx;
            ++spcCtr;
          }
        }
        spcCtr=0;
        for(k=0; k<(int)ckrobj->reactions[rxnIdx].reactants.size(); k++) {

          spcName=ckrobj->reactions[rxnIdx].reactants[k].name;
          spcIdx=ckrobj->speciesData[spcName].index;
          intCount= (int)(ckrobj->reactions[rxnIdx].reactants[k].number+0.5);
	  for(m=0; m<intCount; m++) {
            productSpcIdx[addrProductOfStep[j]+spcCtr]=spcIdx;
            ++spcCtr;
          }
        }  
      }
    } //end if(isNonIntegerStoich[rxnIdx] == 0)

  } // for(j=0; j<nStep; j++)

  return 0;
}

void info_net::setReactionFlags(ckr::CKReader *ckrobj)
{
  int j, k;

  nRxn=static_cast<int>(ckrobj->reactions.size());

  isReversible = new int[nRxn];
  isThirdBody  = new int[nRxn];
  isFalloff    = new int[nRxn]; 
  isNonIntegerStoich = new int[nRxn];

  for(j=0; j<nRxn; j++) {

    isReversible[j]=0;
    if (ckrobj->reactions[j].isReversible && 
	ckrobj->reactions[j].krev.A != 0.0) {
      isReversible[j]=1;
    }

    isThirdBody[j]=0;
    if(ckrobj->reactions[j].isThreeBodyRxn) {
      isThirdBody[j]=1;
    }

    isFalloff[j]=0;
    if(ckrobj->reactions[j].isFalloffRxn) {
      isFalloff[j]=-1;
      if(strcmp(ckrobj->reactions[j].thirdBody.c_str(),"M")!=0) {
        isFalloff[j]=1+spcIdxOfString(*ckrobj,ckrobj->reactions[j].thirdBody);
      }
    }

    const char * all_non_integer_char = getenv("ZERORK_ALL_NON_INTEGER");
    int all_non_integer = 0;
    if (all_non_integer_char != NULL)
    {
       all_non_integer = atoi(all_non_integer_char);
    }
    if( all_non_integer != 0 || ckrobj->reactions[j].isRealOrder) {
      isNonIntegerStoich[j] = 1;
    }
    else {
      isNonIntegerStoich[j] = 0;
      double count;
      int int_count;
      // scan the reactants for non-integer stoichiometric coefficients
      for(k=0; k<(int)ckrobj->reactions[j].reactants.size(); ++k) {
        count = ckrobj->reactions[j].reactants[k].number;
        int_count = (int)(count + 0.5);
        if(fabs(count-(double)int_count) > STOICH_TOL) {
          // set flag that a non-integer stoichiometric coefficient is found
          isNonIntegerStoich[j] = 1;
          break;
        }
      }
      if(isNonIntegerStoich[j] == 0) {
        // scan the products for non-integer stoichiometric coefficients
        for(k=0; k<(int)ckrobj->reactions[j].products.size(); ++k) {
          count = ckrobj->reactions[j].products[k].number;
          int_count = (int)(count + 0.5);
          if(fabs(count-(double)int_count) > STOICH_TOL) {
            // set flag that a non-integer stoichiometric coefficient is found
            isNonIntegerStoich[j] = 1;
            break;
          }
        }
      }
      //
      // TODO - check if reaction is Real Order 
      //
    }
  }   
}


void info_net::print(const species *spcList)
{
  int j,k,spcIdx;
  printf("# From void info_net::print(const species *)\n");
  printf("#\n");
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

int info_net::getSpecIdxOfStepReactant(const int stepId, 
					   const int reacId) const
{
  if(stepId < 0 || stepId >= nStep)
    {return MIN_INT32;}
  if(reacId < 0 || reacId >= nReactantOfStep[stepId])
    {return MIN_INT32;}
  
  return reactantSpcIdx[addrReactantOfStep[stepId]+reacId];
}
int info_net::getSpecIdxOfStepProduct(const int stepId,
					  const int prodId) const
{
  if(stepId < 0 || stepId >= nStep)
    {return MIN_INT32;}
  if(prodId < 0 || prodId >= nProductOfStep[stepId])
    {return MIN_INT32;}

  return productSpcIdx[addrProductOfStep[stepId]+prodId];
}

int info_net::spcIdxOfString(ckr::CKReader &ckrobj, string spcName)
{
  int j;
  int nSpc = static_cast<int>(ckrobj.species.size());
  for(j=0; j<nSpc; j++)
    {
      if(ckrobj.species[j].name == spcName)
	{return j;}
    }
  return zerork::MIN_INT32;
}

} // namespace zerork
