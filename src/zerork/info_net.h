#ifndef ZERORK_INFO_NET_H
#define ZERORK_INFO_NET_H

#include "../CKconverter/CKReader.h"
#include "species.h"
#include "non_integer_reaction_network.h"

namespace zerork {

typedef struct
{
  std::vector<int> index;
  std::vector<double> alpha;

} EnhancedSpecies;

class info_net
{
 public:
  info_net(ckr::CKReader *ckrobj);
  ~info_net();

  int getNumRxn()  const {return nRxn;}
  int getNumSteps() const {return nStep;}
  int getTotalReactants() const {return totReactant;}
  int getTotalProducts() const {return totProduct;}
  int getMaxProductsInStep() const {return maxProductInStep;}
  int getMaxReactantsInStep() const {return maxReactantInStep;}
  int getNumEnhancedSpeciesOfStep(const int id) const 
  {return enhanced_species_of_step[id].index.size();}
  int getEnhancementFactorsOfStep(const int id,
                                  std::vector<int> *species_id,
                                  std::vector<double> *alpha);
  
  int getRxnIdxOfStep(int id) const {return rxnIdxOfStep[id]%nRxn;}
  int getRxnDirOfStep(int id) const
    {return ((rxnIdxOfStep[id]<nRxn) ? 1 : -1);}
  int getStepIdxOfRxn(int id, int dir) const
    {return ((dir > 0) ? stepIdxOfFwdRxn[id] : stepIdxOfRevRxn[id]);}
  // getOrderOfStep and getNumProductsOfStep return zero if the step
  // is a non-integer reaction
  int getOrderOfStep(int id) const {return nReactantOfStep[id];}
  int getNumProductsOfStep(int id) const {return nProductOfStep[id];}
  // getRealOrderOfStep return the real number of the step for exact
  // conversion of the rate constant
  double getRealOrderOfStep(int id) const;
  // bound check safe
  int getSpecIdxOfStepReactant(const int stepId, const int reacId) const;
  int getSpecIdxOfStepProduct(const int stepId, const int prodId) const;

  int getReversibleFlagOfReaction(int id) const {return isReversible[id];} 
  int getThirdBodyFlagOfReaction(int id) const {return isThirdBody[id];} 
  int getFalloffFlagOfReaction(int id) const {return isFalloff[id];} 
  int getNonIntegerStoichFlagOfReaction(int id) const 
    {return isNonIntegerStoich[id];}

  NonIntegerReactionNetwork getNonIntegerReactionNetwork() const 
    {return non_integer_network;}

  void print(const species *spcList);

 private:
  // reaction index <-> step index maps
  int nRxn;
  int *stepIdxOfFwdRxn;
  int *stepIdxOfRevRxn;
  int *rxnIdxOfStep;
  int setRxnStepIdxMaps(ckr::CKReader *ckrobj);
  int SetEnhancedSpeciesList(ckr::CKReader *ckrobj);
  
  // reaction participant species info:
  int nStep;
  int totProduct;
  int totReactant;
  int maxProductInStep;
  int maxReactantInStep;
  int *nReactantOfStep;
  int *nProductOfStep;
  int *addrReactantOfStep;
  int *addrProductOfStep;
  int *reactantSpcIdx;
  int *productSpcIdx;
  
  int *isReversible; // (0,1)
  int *isThirdBody;  // (0,1)
  int *isFalloff;    // (-1 = yes, M, 0 = no, spcId + 1 = yes, falloff species)
  int *isEnhanced;   // (0,1)
  int *isNonIntegerStoich; // (0,1)

  std::vector<EnhancedSpecies> enhanced_species_of_step; 

  NonIntegerReactionNetwork non_integer_network;
 
  int setParticipants(ckr::CKReader *ckrobj);
  int spcIdxOfString(ckr::CKReader &ckrobj, string spcName);
  void setReactionFlags(ckr::CKReader *ckrobj);

 };

} // namespace zerork

#endif
