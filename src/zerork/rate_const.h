#ifndef ZERORK_RATE_CONST_H
#define ZERORK_RATE_CONST_H

#include <string.h>
#include <vector>
#include "../CKconverter/CKReader.h"
#include "info_net.h"
#include "nasa_poly.h"
#include "external_funcs.h"
#include "plog_reaction.h"

namespace zerork {

enum FalloffReactionType {LINDEMANN, 
                          TROE_THREE_PARAMS, 
                          TROE_FOUR_PARAMS,
                          SRI};


typedef struct
{
  int stepIdx;
  double A;
  double Tpow;
  double Tact;
} arrheniusSortElem;

typedef struct
{
  int stepIdx;
  int arrheniusIdx;
} arrheniusStep;

typedef struct
{
  int stepIdx;
  int fwdStepIdx;
  int nReac;
  int nProd;
  double nDelta;
  vector <int> reacSpcIdx;
  vector <int> prodSpcIdx;
} fromKeqStep;

typedef struct
{
  int fwdStepIdx;
  int revStepIdx;
  int nEnhanced;
  vector <int> etbSpcIdx;
  vector <double> etbSpcEff;
} thirdBodyRxn;

typedef struct
{
  int falloffSpcIdx;
  FalloffReactionType falloffType;
  int fwdStepIdx;
  int revStepIdx;
  int nEnhanced;
  vector <int> etbSpcIdx;
  vector <double> etbSpcEff;
  vector <double> param;
} falloffRxn;



int isSameArrheniusTol(arrheniusSortElem x, arrheniusSortElem y);
int compareArrhenius(const void *x, const void *y); 
int compareArrheniusT1000(const void *x, const void *y);

class rate_const
{
 public:
  rate_const(ckr::CKReader *ckrobj, info_net *netobj, nasa_poly_group *tobj);
  virtual ~rate_const();

  void updateK(const double T, const double C[]);
  void updateK(const double T, const double C[], double Kcopy[]);
  void updateK_TCM(const double T, const double C[], const double C_mix);
  void updateK_TCM(const double T, const double C[], 
                   const double C_mix, double Kcopy[]);
  void updateKExplicit(const double T, const double C[], double Kcopy[]);
  double * getKptr() const {return &Kwork[0];}
  void getKrxn(info_net &netobj, double Kfwd[], double Krev[]);
  int getNumSpecies() const {return nSpc;}
  void print();

//  void writeExplicitUpdates(const char *, const char *);
  void write_funcs(FILE* fptr);
  void setUseExArrh(){ use_external_arrh = true; };
  void unsetUseExArrh(){ use_external_arrh = false; };
  void setUseExKeq(){ use_external_keq = true; };
  void unsetUseExKeq(){ use_external_keq = false; };
  void setExArrhFunc(external_func_arrh_t fn_handle) { ex_func_calc_arrh = fn_handle; };
  void setExKeqFunc(external_func_keq_t fn_handle) { ex_func_calc_keq = fn_handle; };

 protected:
 
  int nSpc;
  int nStep;
  int cpySize;
  double *Kwork;    // length nStep
  double *Gibbs_RT; // length nSpc 

  double convertE;
  double convertC;

  double Csum;
  bool Tchanged;
  double Tcurrent;
  double log_e_Tcurrent;
  double invTcurrent;
  double log_e_PatmInvRuT;
  void updateTcurrent(double const T);

  // sizes of temperature based terms
  int nArrheniusStep;   
  int nLandauTellerStep;
  int nFromKeqStep;
  int nPLogInterpolationStep;
  void setStepCount_Ttype(ckr::CKReader *ckrobj);
  // nStep must equal nArrheniusStep+nLandauTellerStep+nFromKeqStep

  // sizes of pressure based terms 
  int nFalloffRxn;
  int nThirdBodyRxn;
  void setRxnCount_Ptype(ckr::CKReader &ckrobj);
  thirdBodyRxn *thirdBodyRxnList;
  void setThirdBodyRxnList(ckr::CKReader &ckrobj, info_net &netobj);
  void updateThirdBodyRxn(const double C[]);
  int spcIdxOfString(ckr::CKReader &ckrobj, string spcName);
  int getThirdBodyEff(ckr::CKReader &ckrobj, int rxnId, vector <int> &spcId,
		      vector <double> &spcEff);
  falloffRxn *falloffRxnList;
  void setFalloffRxnList(ckr::CKReader &ckrobj, info_net &netobj);
  void updateFalloffRxn(const double C[]);
  int isNonStandardTroe(const int falloffId, const int rxnId) const;

  std::vector<PLogReaction> plogInterpolationStepList; 
  void setPLogInterpolationStepList(ckr::CKReader &ckrobj, info_net &netobj);
  void updatePLogInterpolationStep(const double pressure, 
                                   const double log_e_pressure);

  // 
  int nDistinctArrhenius;
  arrheniusStep *arrheniusStepList;
  double *distinctArrheniusLogAfact;
  double *distinctArrheniusTpow;
  double *distinctArrheniusTact;
  double *arrheniusCoeffs;
  double *arrWorkArray;
  double *keqWorkArray;
  double *falloffWorkArray;
  void setArrheniusStepList(ckr::CKReader *ckrobj, info_net *netobj);
  void updateArrheniusStep();

  fromKeqStep *fromKeqStepList;
  void setFromKeqStepList(ckr::CKReader &ckrobj, info_net &netobj);
  void updateFromKeqStep();

  nasa_poly_group *thermoPtr;

  bool use_external_arrh;
  bool use_external_keq;
  external_func_arrh_t ex_func_calc_arrh;
  external_func_keq_t ex_func_calc_keq;

  bool use_non_integer_network_;
  NonIntegerReactionNetwork non_integer_network_;
};

} // namespace zerork

#endif
