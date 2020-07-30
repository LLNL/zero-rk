#ifndef ZERORK_SIMPLE_NETWORK_H
#define ZERORK_SIMPLE_NETWORK_H


#include "../CKconverter/CKReader.h"
#include "nasa_poly.h"
#include "species.h"

namespace zerork {

class simple_network
{
 public:
  simple_network(ckr::CKReader *ckrobj, nasa_poly_group &thermo);
  ~simple_network();

  void print_info(const species *spcList);

 private:
  double convertC;
  double convertE; 

  // reaction index <-> step index maps
  int nRxn;
  int *stepIdxOfFwdRxn;
  int *stepIdxOfRevRxn;
  int *rxnIdxOfStep;
  int setRxnStepIdxMaps(ckr::CKReader *ckrobj);
    
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
  int setParticipants(ckr::CKReader *ckrobj);
 };

} // namespace zerork

#endif
