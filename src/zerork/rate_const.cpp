#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <algorithm> // for std::max

#include "rate_const.h"
#include "constants.h"
#include "fast_exps.h"

#ifndef HAVE_ALIGNED_ALLOC
#if defined(__APPLE__)
static inline void* aligned_alloc(size_t alignment, size_t size)
{
        void* p;
        int flag = posix_memalign(&p, alignment, size);
        return (flag == 0) ? p : 0;
}
#elif defined(_WIN32)
static inline void* aligned_alloc(size_t alignment, size_t size)
{
        void* p = _aligned_malloc(size, alignment);
        return (p == nullptr) ? 0 : p;
}
#endif
#endif
#ifndef _WIN32
#define _aligned_free free
#endif

#ifndef _WIN32
#define unlikely(expr) __builtin_expect(!!(expr), 0)
#define likely(expr) __builtin_expect(!!(expr), 1)
#else
#define unlikely(expr) expr
#define likely(expr) expr
#endif

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

rate_const::rate_const(ckr::CKReader *ckrobj, info_net *netobj,
			       nasa_poly_group *tobj)
{
  assert(ckrobj != NULL);
  assert(netobj != NULL);
  //N.B. all checks for validity should be done in CKReader/CKParser.
  //  asserts are used here to make sure that rate_const is kept up
  //  with CKReader, but they should not be catching issues with input
  //  files, such that if we get here we have parsed the inputs and
  //  know that they are OK.

  // determine the multiplier to convert the activation energy units in the
  // prescribed mech into kelvin
  convertE = 0;
  if(ckrobj->units.ActEnergy == ckr::Cal_per_Mole)
    {convertE = CAL_PER_MOL_TACT;}
  else if(ckrobj->units.ActEnergy == ckr::Kcal_per_Mole)
    {convertE = KCAL_PER_MOL_TACT;}
  else if(ckrobj->units.ActEnergy == ckr::Kelvin)
    {convertE = 1.0;}
  else if(ckrobj->units.ActEnergy == ckr::Kjoules_per_Mole)
    {convertE = KJOULES_PER_MOL_TACT;}
  else if(ckrobj->units.ActEnergy == ckr::Joules_per_Mole)
    {convertE = JOULES_PER_MOL_TACT;}
  assert(convertE!=0); //This is to make sure we handled parsing correctly in ckr

  // convert the concentration units in the rate
  // constants to kmol/m^3
  assert(ckrobj->units.Quantity == ckr::Moles);
  // note that the length scale in the CKReader object is [cm]
  // [mol/cm^3] * convertC = [kmol/m^3]
  convertC=1000.0;

  thermoPtr=tobj;
  nStep=netobj->getNumSteps();
  cpySize=nStep*sizeof(double);
  Kwork = new double[nStep];
  nSpc=ckrobj->species.size();
  Gibbs_RT = new double[nSpc];

  setStepCount_Ttype(ckrobj); // count tempertature dependent step types
  setRxnCount_Ptype(*ckrobj); // count pressure     dependent rxn  types
  setArrheniusStepList(ckrobj,netobj);

  // non_integer network must be set up before the rate constants from
  // chemical equilibrium parameters are set up
  non_integer_network_ = netobj->getNonIntegerReactionNetwork();
  use_non_integer_network_ = true;
  if(non_integer_network_.GetNumNonIntegerSteps() == 0) {
    use_non_integer_network_ = false;
  }
  setFromKeqStepList(*ckrobj,*netobj);
  setThirdBodyRxnList(*ckrobj,*netobj);
  setFalloffRxnList(*ckrobj,*netobj);
  setPLogInterpolationStepList(*ckrobj,*netobj);

  Tchanged = true;
  Tcurrent = 0;

  int allocSize = nDistinctArrhenius*sizeof(double);
  allocSize = ((allocSize + 31)/32)*32; //round to next even multiple of 4 doubles
  arrWorkArray = nullptr;
  if(allocSize > 0) {
    arrWorkArray = (double*)aligned_alloc(32, allocSize);
    memset(arrWorkArray,0.0, allocSize);
  }
  allocSize = nFromKeqStep*sizeof(double);
  allocSize = ((allocSize + 31)/32)*32; //round to next even multiple of 4 doubles
  keqWorkArray = nullptr;
  if(allocSize > 0) {
    keqWorkArray = (double*)aligned_alloc(32, allocSize);
    memset(keqWorkArray, 0.0, allocSize);
  }
  allocSize = nFalloffRxn*sizeof(double);
  allocSize = ((allocSize + 31)/32)*32; //round to next even multiple of 4 doubles
  falloffWorkArray = nullptr;
  if(allocSize > 0) {
    falloffWorkArray = (double*)aligned_alloc(32, allocSize);
    memset(falloffWorkArray, 0.0, allocSize);
  }

  use_external_arrh = false;
  use_external_keq = false;

}

rate_const::~rate_const()
{
  delete [] falloffRxnList;
  delete [] thirdBodyRxnList;
  delete [] fromKeqStepList;
  if(arrheniusStepList != NULL) {
    delete [] arrheniusStepList;
    delete [] distinctArrheniusLogAfact;
    delete [] distinctArrheniusTpow;
    delete [] distinctArrheniusTact;
  }
  delete [] Gibbs_RT;
  delete [] Kwork;
  if(arrWorkArray != nullptr) {
    _aligned_free(arrWorkArray);
  }
  if(keqWorkArray != nullptr) {
    _aligned_free(keqWorkArray);
  }
  if(falloffWorkArray != nullptr) {
    _aligned_free(falloffWorkArray);
  }
}

void rate_const::setStepCount_Ttype(ckr::CKReader *ckrobj)
{
  int j;
  int nRxn=static_cast<int>(ckrobj->reactions.size());
  nArrheniusStep=0;
  nLandauTellerStep=0;
  nFromKeqStep=0;
  nPLogInterpolationStep=0;

  for(j=0; j<nRxn; j++)
    {
      if(ckrobj->reactions[j].kf.type == ckr::Arrhenius)
	{++nArrheniusStep;}
      else if(ckrobj->reactions[j].kf.type == ckr::PLogInterpolation)
        {++nPLogInterpolationStep;}
      else {
        assert(("Unsupported forward reaction type", false));
      }

      if(ckrobj->reactions[j].isReversible && ckrobj->reactions[j].krev.A != 0.0)
	{
          // krev.type =ckr::Arrhenius for REV keyword and for reversible
          // reactions with no reverse parameters specified, including PLOG
	  assert(("Unsupported reverse reaction type.", ckrobj->reactions[j].krev.type == ckr::Arrhenius));
	  if(ckrobj->reactions[j].krev.A > 0.0)
	    {++nArrheniusStep;}
	  else
	    {++nFromKeqStep;}
 	}
    }
}

void rate_const::setRxnCount_Ptype(ckr::CKReader &ckrobj)
{
  int j;
  int nRxn=static_cast<int>(ckrobj.reactions.size());
  nFalloffRxn=nThirdBodyRxn=0;
  for(j=0; j<nRxn; j++) {
    if(ckrobj.reactions[j].isFalloffRxn) {
      ++nFalloffRxn;
    }
    if(ckrobj.reactions[j].isThreeBodyRxn) {
      ++nThirdBodyRxn;
    }
    if(ckrobj.reactions[j].isFalloffRxn &&
       ckrobj.reactions[j].isThreeBodyRxn) {
      assert(("Can't be both fall-off and 3rd-body", false));
    }
  }
}


void rate_const::setArrheniusStepList(ckr::CKReader *ckrobj,
					  info_net *netobj)
{
  int j,k,rxnIdx,rxnDir;
  arrheniusSortElem *paramListTmp;
  arrheniusSortElem lastDistinct;

  if(nArrheniusStep > 0) {

    paramListTmp = new arrheniusSortElem[nArrheniusStep];
    k=0; // Arrhenius step counter

    for(j=0; j<nStep; j++) {
      rxnIdx=netobj->getRxnIdxOfStep(j);
      rxnDir=netobj->getRxnDirOfStep(j);

      if(rxnDir == 1) { // forward
        if(ckrobj->reactions[rxnIdx].kf.type == ckr::Arrhenius) {
	  paramListTmp[k].stepIdx=j;
          // convert and store the arrhenius parameters
	  paramListTmp[k].A = ckrobj->reactions[rxnIdx].kf.A*
	    pow(convertC,(1.0-netobj->getRealOrderOfStep(j)));

	  paramListTmp[k].Tpow = ckrobj->reactions[rxnIdx].kf.n;

	  paramListTmp[k].Tact = ckrobj->reactions[rxnIdx].kf.E*convertE;

	  // if the reaction is a third body/non-falloff type reaction
	  // an additional concentration factor is needed for the
	  // thirdbody term
	  if(ckrobj->reactions[rxnIdx].isThreeBodyRxn)
	    {paramListTmp[k].A/=convertC;}
          ++k; // increment Arrhenius counter
        }
      }
      else { // reverse direction
	if(ckrobj->reactions[rxnIdx].krev.type == ckr::Arrhenius &&
           ckrobj->reactions[rxnIdx].krev.A > 0.0)
	  {
	    paramListTmp[k].stepIdx=j;
	    // convert and store the arrhenius parameters
	    paramListTmp[k].A = ckrobj->reactions[rxnIdx].krev.A*
	      pow(convertC,(1.0-netobj->getRealOrderOfStep(j)));

	    paramListTmp[k].Tpow = ckrobj->reactions[rxnIdx].krev.n;

	    paramListTmp[k].Tact = ckrobj->reactions[rxnIdx].krev.E*convertE;

	    // if the reaction is a third body/non-falloff type reaction
	    // an additional concentration factor is needed for the
	    // thirdbody term
	    if(ckrobj->reactions[rxnIdx].isThreeBodyRxn)
	      {paramListTmp[k].A/=convertC;}
            ++k; // increment Arrhenius counter
	  }
       }
    } // end for(j=0; j<nStep; j++)

    //The number of arrheniusSteps recorded (k)
    //should match the same as the number found in
    //rate_const::setStepCount_Ttype(...) (nArrheniusStep)
    assert(("Logical error in step counts.", k == nArrheniusStep));

    qsort((void *)paramListTmp,nArrheniusStep,sizeof(arrheniusSortElem),
          compareArrhenius);
    //for(j=0; j<nStep; j++)
    //  {
    //    printf("rank %d: step = %d (A,Tpow,Tact) = (%16.8e,%8.5f,%10.2f)\n",
    //	     j,paramListTmp[j].stepIdx,paramListTmp[j].A,paramListTmp[j].Tpow,
    //	     paramListTmp[j].Tact);
    //}

    nDistinctArrhenius=1;
    lastDistinct=paramListTmp[0];
    for(j=1; j<nArrheniusStep; j++) {

      if(isSameArrheniusTol(lastDistinct,paramListTmp[j])==0) {
	lastDistinct=paramListTmp[j];
	++nDistinctArrhenius;
      }
    }

    distinctArrheniusLogAfact = new double[nDistinctArrhenius];
    distinctArrheniusTpow  = new double[nDistinctArrhenius];
    distinctArrheniusTact  = new double[nDistinctArrhenius];
    arrheniusStepList = new arrheniusStep[nArrheniusStep];

    distinctArrheniusLogAfact[0]=log(paramListTmp[0].A);
    distinctArrheniusTpow[0]=paramListTmp[0].Tpow;
    distinctArrheniusTact[0]=paramListTmp[0].Tact;

    arrheniusStepList[0].arrheniusIdx=0;
    arrheniusStepList[0].stepIdx=paramListTmp[0].stepIdx;

    k=1; // distinct arrhenius counter
    lastDistinct=paramListTmp[0];
    for(j=1; j<nArrheniusStep; j++)
    {
      // new arrhenius expression found
      if(isSameArrheniusTol(lastDistinct,paramListTmp[j])==0)
	{
	  distinctArrheniusLogAfact[k]=log(paramListTmp[j].A);
	  distinctArrheniusTpow[k]=paramListTmp[j].Tpow;
	  distinctArrheniusTact[k]=paramListTmp[j].Tact;

	  lastDistinct=paramListTmp[j]; ++k;
	}

      arrheniusStepList[j].arrheniusIdx=k-1;
      arrheniusStepList[j].stepIdx=paramListTmp[j].stepIdx;
    }
    delete [] paramListTmp;
  } // end if(nArrheniusStep > 0)
  else {
    nDistinctArrhenius=0;
    distinctArrheniusLogAfact = NULL;
    distinctArrheniusTpow     = NULL;
    distinctArrheniusTact     = NULL;
    arrheniusStepList         = NULL;
  }
}

void rate_const::setFromKeqStepList(ckr::CKReader &ckrobj,
					info_net &netobj)
{
  int j,k,m,fwdIdx,revIdx;
  int nRxn=netobj.getNumRxn();

  fromKeqStepList = new fromKeqStep[nFromKeqStep];

  k=0; // fromKeq counter
  for(j=0; j<nRxn; j++)
    {
      if(ckrobj.reactions[j].isReversible && ckrobj.reactions[j].krev.A < 0.0)
	{
	  fwdIdx=netobj.getStepIdxOfRxn(j,1);
	  revIdx=netobj.getStepIdxOfRxn(j,-1);

	  fromKeqStepList[k].stepIdx=revIdx;
	  fromKeqStepList[k].fwdStepIdx=fwdIdx;
          // If reaction j is a non-integer reaction, then nReac and nProd
          // are zero.
	  fromKeqStepList[k].nReac=netobj.getOrderOfStep(fwdIdx);
	  fromKeqStepList[k].nProd=netobj.getNumProductsOfStep(fwdIdx);
	  fromKeqStepList[k].nDelta=(double)(fromKeqStepList[k].nProd-
					     fromKeqStepList[k].nReac);
          if(use_non_integer_network_) {
            if(non_integer_network_.HasReaction(j)) {
              // number of products as defined by the reverse rate-of-progress
              // concentration powers minus the number of reactants as
              // defined by the forward rate-of-progress concentration powers
              fromKeqStepList[k].nDelta = // nProd - nReac
                non_integer_network_.GetOrderOfStep(revIdx) -
                non_integer_network_.GetOrderOfStep(fwdIdx);
            }
          }
          // If reaction j is a non-integer reaction, then nReac and nProd
          // are zero and the following assignments are skipped
	  fromKeqStepList[k].reacSpcIdx.resize(fromKeqStepList[k].nReac);
	  fromKeqStepList[k].prodSpcIdx.resize(fromKeqStepList[k].nProd);

	  for(m=0; m<fromKeqStepList[k].nReac; m++)
	    {
	      fromKeqStepList[k].reacSpcIdx[m]=
		netobj.getSpecIdxOfStepReactant(fwdIdx,m);
	    }
	  for(m=0; m<fromKeqStepList[k].nProd; m++)
	    {
	      fromKeqStepList[k].prodSpcIdx[m]=
		netobj.getSpecIdxOfStepProduct(fwdIdx,m);
	    }
	  ++k;
	}
    }
}

void rate_const::setThirdBodyRxnList(ckr::CKReader &ckrobj,
					 info_net &netobj)
{
  int j,k,nEnh;
  int nRxn=netobj.getNumRxn();
  string spcName;

  thirdBodyRxnList=new thirdBodyRxn[nThirdBodyRxn];

  k=0; // thirdBodyRxn counter
  for(j=0; j<nRxn; j++) {

    if(ckrobj.reactions[j].isThreeBodyRxn) {

      thirdBodyRxnList[k].fwdStepIdx=netobj.getStepIdxOfRxn(j,1);
      thirdBodyRxnList[k].revStepIdx=netobj.getStepIdxOfRxn(j,-1);
      nEnh=ckrobj.reactions[j].e3b.size();

      if(nEnh > 0) {

        thirdBodyRxnList[k].etbSpcIdx.resize(nEnh);
	thirdBodyRxnList[k].etbSpcEff.resize(nEnh);
	nEnh=getThirdBodyEff(ckrobj,j,thirdBodyRxnList[k].etbSpcIdx,
			     thirdBodyRxnList[k].etbSpcEff);
	thirdBodyRxnList[k].etbSpcIdx.resize(nEnh);
	thirdBodyRxnList[k].etbSpcEff.resize(nEnh);
      }

      thirdBodyRxnList[k].nEnhanced=nEnh;
      ++k;
    }
  }
}

void rate_const::setFalloffRxnList(ckr::CKReader &ckrobj,
				       info_net &netobj)
{
  int j,k,nEnh;
  int nRxn=netobj.getNumRxn();
  double log_e_Alow;
  string spcName;

  falloffRxnList=new falloffRxn[nFalloffRxn];

  k=0; // falloffRxn counter
  for(j=0; j<nRxn; j++) {

    if(ckrobj.reactions[j].isFalloffRxn) {

      falloffRxnList[k].fwdStepIdx=netobj.getStepIdxOfRxn(j,1);
      falloffRxnList[k].revStepIdx=netobj.getStepIdxOfRxn(j,-1);
      nEnh=ckrobj.reactions[j].e3b.size();

      falloffRxnList[k].falloffSpcIdx=MIN_INT32;
      if(strcmp(ckrobj.reactions[j].thirdBody.c_str(),"M")!=0) {

        falloffRxnList[k].falloffSpcIdx=
          spcIdxOfString(ckrobj,ckrobj.reactions[j].thirdBody);
      }

      if(nEnh > 0) {

        falloffRxnList[k].etbSpcIdx.resize(nEnh);
	falloffRxnList[k].etbSpcEff.resize(nEnh);
	nEnh=getThirdBodyEff(ckrobj,j,falloffRxnList[k].etbSpcIdx,
			     falloffRxnList[k].etbSpcEff);
	falloffRxnList[k].etbSpcIdx.resize(nEnh);
	falloffRxnList[k].etbSpcEff.resize(nEnh);
      }

      if(nEnh > 0 &&
         falloffRxnList[k].falloffSpcIdx >= 0) {
        // enhanced third body reactions ignored for single falloff species
        printf("# INFO: reaction %d has a single falloff species (+%s)\n",
               j,ckrobj.reactions[j].thirdBody.c_str());
        printf("#       species index %d, and %d enhanced third body\n",
               falloffRxnList[k].falloffSpcIdx,nEnh);
        printf("#       species enhancements defined.  The third body\n");
        printf("#       enhancements will be ignored.\n");
      }

      falloffRxnList[k].nEnhanced=nEnh;

      falloffRxnList[k].param.resize(7);
      log_e_Alow = ckrobj.reactions[j].kf_aux.A*
        pow(convertC,
	   -netobj.getRealOrderOfStep(falloffRxnList[k].fwdStepIdx));

      falloffRxnList[k].param[0]=log(log_e_Alow);
      falloffRxnList[k].param[1]=ckrobj.reactions[j].kf_aux.n;
      falloffRxnList[k].param[2]=ckrobj.reactions[j].kf_aux.E*convertE;

      if(ckrobj.reactions[j].falloffType == ckr::Lindemann) {

        falloffRxnList[k].falloffType = LINDEMANN;
        falloffRxnList[k].param.resize(3);
	//printf(" *** rxn %d Lindemann falloff\n",j+1);
      }
      else if(ckrobj.reactions[j].falloffType == ckr::Troe) {
        // alpha
        falloffRxnList[k].param[3]=ckrobj.reactions[j].falloffParameters[0];
	// T***
	falloffRxnList[k].param[4]=ckrobj.reactions[j].falloffParameters[1];
	// T*
        falloffRxnList[k].param[5]=ckrobj.reactions[j].falloffParameters[2];

        int nParams = ckrobj.reactions[j].falloffParameters.size();
        assert(("Incorrect number of falloff parameters", nParams==3 || nParams==4));
        if(nParams==3) {
          // 3-parameter TROE
          falloffRxnList[k].falloffType = TROE_THREE_PARAMS;
          falloffRxnList[k].param.resize(6);
        }
        else if(nParams==4) {
          // 4-parameter TROE
          falloffRxnList[k].falloffType = TROE_FOUR_PARAMS;
          // T**
          falloffRxnList[k].param[6]=ckrobj.reactions[j].falloffParameters[3];
        }
      } // end of if reaction j is a 3 or 4 parameter Troe reaction
      else if(ckrobj.reactions[j].falloffType == ckr::SRI) {
        //printf("# DEBUG: reaction %d identified as SRI falloff\n",
        //       j);
	//printf("#        with %d parameters\n",
        //       (int)ckrobj.reactions[j].falloffParameters.size());

        falloffRxnList[k].falloffType = SRI;
        falloffRxnList[k].param.resize(3+ckrobj.reactions[j].falloffParameters.size());
        falloffRxnList[k].param[3] = ckrobj.reactions[j].falloffParameters[0];
        falloffRxnList[k].param[4] = ckrobj.reactions[j].falloffParameters[1];
        if(ckrobj.reactions[j].falloffParameters[2] > 0) {
            falloffRxnList[k].param[5] = 1.0/ckrobj.reactions[j].falloffParameters[2];
        } else {
            falloffRxnList[k].param[5] = -1.0;
        }

        // copy the additional parameters if present
        for(size_t m=3; m<ckrobj.reactions[j].falloffParameters.size(); ++m) {
          falloffRxnList[k].param[3+m] =
            ckrobj.reactions[j].falloffParameters[m];
        }
      }
      else {
        assert(("Unsupported falloff reaction type", false));
      }

      ++k;  // falloffRxn counter
    }
  }
}

void rate_const::updateK(const double T, const double C[])
{
  int j;
  double pressure;
  // initialize to aid in debugging
  for(j=0; j<nStep; j++)
    {Kwork[j]=0.0;}

  updateTcurrent(T);
  Csum=0.0;
  for(j=0; j<nSpc;)
    {Csum+=C[j]; ++j;}
  pressure = Csum*NIST_RU*T;

  if(use_external_arrh)
  {
     ex_func_calc_arrh(Tcurrent,arrWorkArray,Kwork,
                       nDistinctArrhenius,
                       distinctArrheniusLogAfact,
                       distinctArrheniusTpow,
                       distinctArrheniusTact);
  }
  else
  {
      updateArrheniusStep();
  }
  // PLOG reactions must be updated before computing the reverse rates from
  // Keq
  updatePLogInterpolationStep(pressure,
                              log(pressure));

  if(use_external_keq)
  {
     thermoPtr->getG_RT(Tcurrent,Gibbs_RT);
     ex_func_calc_keq(nFromKeqStep,Gibbs_RT,keqWorkArray,Kwork,log_e_PatmInvRuT);
  }
  else
  {
      updateFromKeqStep();
  }

  updateThirdBodyRxn(&C[0]);
  updateFalloffRxn(&C[0]);
}

// Update the reaction rate constants using the TCM state variable
// (C = concentration, M = mixture concentration, T = temperature).  The
// function is a copied from updateK, with the Csum calcuation replaced with
// the function argument C_mix.  This allows for the mixture concentration,
// and, in effect, the pressure to be changed independently of the species
// composition concentration.
void rate_const::updateK_TCM(const double T,
                             const double C[],
                             const double C_mix)
{
  int j;
  double pressure;
  // initialize to aid in debugging
  for(j=0; j<nStep; j++)
    {Kwork[j]=0.0;}

  updateTcurrent(T);
  // Original from updateK takes the concetration sum from all the species
  //Csum=0.0;
  //for(j=0; j<nSpc;)
  //  {Csum+=C[j]; ++j;}
  //
  // In updateK_CMT, the C_mix argument is used instead.
  Csum = C_mix;
  pressure = Csum*NIST_RU*T;


  if(use_external_arrh)
  {
     ex_func_calc_arrh(Tcurrent,arrWorkArray,Kwork,
                       nDistinctArrhenius,
                       distinctArrheniusLogAfact,
                       distinctArrheniusTpow,
                       distinctArrheniusTact);
  }
  else
  {
     updateArrheniusStep();
  }
  // PLOG reactions must be updated before computing the reverse rates from
  // Keq
  updatePLogInterpolationStep(pressure,
                              log(pressure));

  if(use_external_keq)
  {
     thermoPtr->getG_RT(Tcurrent,Gibbs_RT);
     ex_func_calc_keq(nFromKeqStep,Gibbs_RT,keqWorkArray,Kwork,log_e_PatmInvRuT);
  }
  else
  {
      updateFromKeqStep();
  }

  updateThirdBodyRxn(&C[0]);
  updateFalloffRxn(&C[0]);
}


void rate_const::updateK(const double T, const double C[], double Kcopy[])
{
  updateK(T,C);
  memcpy(Kcopy,Kwork,cpySize);
}
void rate_const::updateK_TCM(const double T,
                             const double C[],
                             const double C_mix,
                             double Kcopy[])
{
  updateK_TCM(T,C,C_mix);
  memcpy(Kcopy,Kwork,cpySize);
}
void rate_const::updateTcurrent(double const T)
{
  Tchanged=false;
  if(T!=Tcurrent) {
    Tchanged = true;
    Tcurrent=T;
    log_e_Tcurrent=log(Tcurrent);
    invTcurrent=1.0/Tcurrent;
    log_e_PatmInvRuT=log(P_ATM/(NIST_RU*T));
  }
}


void rate_const::updateArrheniusStep()
{
  int j;
  if(Tchanged) {
    //Need below def's for gcc to vectorize the loop
    const double local_log_e_Tcurrent = log_e_Tcurrent;
    const double local_invTcurrent = invTcurrent;
    for(j=0; j<nDistinctArrhenius; ++j) {
        arrWorkArray[j]=distinctArrheniusLogAfact[j]
                 	     +distinctArrheniusTpow[j]*local_log_e_Tcurrent
     	             -distinctArrheniusTact[j]*local_invTcurrent;
    }
    fast_vec_exp(arrWorkArray,nDistinctArrhenius);
  }
  for(j=0; j<nArrheniusStep; ++j) {
      Kwork[arrheniusStepList[j].stepIdx] =
	arrWorkArray[arrheniusStepList[j].arrheniusIdx];
    }
}


void rate_const::updateFromKeqStep()
{
  int j,k;
  if(Tchanged) {
    double thermo_sum=0.0;

    thermoPtr->getG_RT(Tcurrent,Gibbs_RT);

    for(j=0; j<nFromKeqStep; ++j) {

      const int forward_step_id = fromKeqStepList[j].fwdStepIdx;

      if(non_integer_network_.HasStep(forward_step_id)) {
        // products - reactants (defined relative to the forward direction)
        thermo_sum =
          non_integer_network_.GetThermoChangeOfStep(forward_step_id,
                                                     Gibbs_RT);

      } else {
        // the reactant and product counts are defined relative to the forward
        // step direction
        const int num_reactants = fromKeqStepList[j].nReac;
        const int num_products  = fromKeqStepList[j].nProd;
        thermo_sum=0.0;
        for(k=0; k<num_products; ++k) {
          thermo_sum += Gibbs_RT[fromKeqStepList[j].prodSpcIdx[k]];
        }
        for(k=0; k<num_reactants; ++k) {
          thermo_sum -= Gibbs_RT[fromKeqStepList[j].reacSpcIdx[k]];
        }

      }
      keqWorkArray[j] = thermo_sum-fromKeqStepList[j].nDelta*log_e_PatmInvRuT;
    }
    fast_vec_exp(keqWorkArray,nFromKeqStep);
  }
  for(j=0; j<nFromKeqStep; j++) {

    Kwork[fromKeqStepList[j].stepIdx]=
      keqWorkArray[j]*Kwork[fromKeqStepList[j].fwdStepIdx];
  }
}

//void rate_const::writeExplicitUpdates(const char *fileName, const char *funcName)
void rate_const::write_funcs(FILE *fptr)
{
  UnsupportedFeature(__FILE__,__LINE__);
  int j,k;

  fprintf(fptr,"void external_func_arrh(const double Tcurrent, double arrWorkArray[], double Kwork[],\n");
  fprintf(fptr,"                            const int nDistinctArrhenius,\n");
  fprintf(fptr,"                            const double distinctArrheniusLogAfact[],\n");
  fprintf(fptr,"                            const double distinctArrheniusTpow[],\n");
  fprintf(fptr,"                            const double distinctArrheniusTact[])\n");
  fprintf(fptr,"{\n");
#ifdef ZERORK_ARRH_FUNC_LOOP
  fprintf(fptr,"    int j;\n\n");
  fprintf(fptr,"    const double local_log_e_Tcurrent = log(Tcurrent);\n");
  fprintf(fptr,"    const double local_invTcurrent = 1.0/Tcurrent;\n");
  fprintf(fptr,"    for(j=0; j<nDistinctArrhenius; ++j)\n");
  fprintf(fptr,"      {\n");
  fprintf(fptr,"        arrWorkArray[j]= distinctArrheniusLogAfact[j]\n");
  fprintf(fptr,"                        +distinctArrheniusTpow[j]*local_log_e_Tcurrent\n");
  fprintf(fptr,"                        -distinctArrheniusTact[j]*local_invTcurrent;\n");
  fprintf(fptr,"      }\n");
#else
  fprintf(fptr,"    int j;\n\n");
  fprintf(fptr,"    const double local_log_e_Tcurrent = log(Tcurrent);\n");
  fprintf(fptr,"    const double local_invTcurrent = 1.0/Tcurrent;\n");
  fprintf(fptr,"\n");

  fprintf(fptr,"    static const double LogAfact[%d] = {%a",nDistinctArrhenius,distinctArrheniusLogAfact[0]);
  for(j=1; j<nDistinctArrhenius; ++j)
  {
      if((j+1) % 5 == 0)
      {
          fprintf(fptr,",\n                                          ");
      }
      else
      {
          fprintf(fptr,", ");
      }
      fprintf(fptr,"%a",distinctArrheniusLogAfact[j]);
  }
  fprintf(fptr,"};\n");
  fprintf(fptr,"\n");

  fprintf(fptr,"    static const double Tpow[%d] = {%a",nDistinctArrhenius,distinctArrheniusTpow[0]);
  for(j=1; j<nDistinctArrhenius; ++j)
  {
      if((j+1) % 5 == 0)
      {
          fprintf(fptr,",\n                                          ");
      }
      else
      {
          fprintf(fptr,", ");
      }
      fprintf(fptr,"%a",distinctArrheniusTpow[j]);
  }
  fprintf(fptr,"};\n");
  fprintf(fptr,"\n");

  fprintf(fptr,"    static const double Tact[%d] = {%a",nDistinctArrhenius,distinctArrheniusTact[0]);
  for(j=1; j<nDistinctArrhenius; ++j)
  {
      if((j+1) % 5 == 0)
      {
          fprintf(fptr,",\n                                          ");
      }
      else
      {
          fprintf(fptr,", ");
      }
      fprintf(fptr,"%a",distinctArrheniusTact[j]);
  }
  fprintf(fptr,"};\n");
  fprintf(fptr,"\n");

  fprintf(fptr,"    for(j=0; j<%d; ++j)\n",nDistinctArrhenius);
  fprintf(fptr,"      {\n");
  fprintf(fptr,"        arrWorkArray[j]= LogAfact[j]\n");
  fprintf(fptr,"                        +Tpow[j]*local_log_e_Tcurrent\n");
  fprintf(fptr,"                        -Tact[j]*local_invTcurrent;\n");
  fprintf(fptr,"      }\n");
#endif


  fprintf(fptr,"\n");
  fprintf(fptr,"    fast_vec_exp(arrWorkArray,%d);\n",(nDistinctArrhenius+nDistinctArrhenius%4));
  fprintf(fptr,"\n");
  for(j=0; j<nArrheniusStep; ++j)
    {
      fprintf(fptr,"    Kwork[%5d]=arrWorkArray[%5d];\n",arrheniusStepList[j].stepIdx,
                                                       arrheniusStepList[j].arrheniusIdx);
    }
  fprintf(fptr,"}\n");
  fprintf(fptr,"\n\n\n");

  fprintf(fptr,"void external_func_keq(const int nFromKeqStep, const double Gibbs_RT[], double keqWorkArray[], double Kwork[], const double log_e_PatmInvRuT)\n");
  fprintf(fptr,"{\n");

  for(j=0; j<nFromKeqStep; j++)
    {
      fprintf(fptr,"    keqWorkArray[%5d]=-( (",j);
      fprintf(fptr,"Gibbs_RT[%5d]",fromKeqStepList[j].reacSpcIdx[0]);
      for(k=1; k<fromKeqStepList[j].nReac; ++k)
	{
          fprintf(fptr,"+Gibbs_RT[%5d]",fromKeqStepList[j].reacSpcIdx[k]);
        }
      fprintf(fptr,") - (");
      fprintf(fptr,"Gibbs_RT[%5d]",fromKeqStepList[j].prodSpcIdx[0]);
      for(k=1; k<fromKeqStepList[j].nProd; ++k)
	{
          fprintf(fptr,"+Gibbs_RT[%5d]",fromKeqStepList[j].prodSpcIdx[k]);
        }
      fprintf(fptr,") + %18.10g*log_e_PatmInvRuT);\n",fromKeqStepList[j].nDelta);
    }
  fprintf(fptr,"\n");
  fprintf(fptr,"    fast_vec_exp(keqWorkArray,%d);\n",(nFromKeqStep+nFromKeqStep%4));
  fprintf(fptr,"\n");
  for(j=0; j<nFromKeqStep; j++)
    {
      fprintf(fptr,"    Kwork[%5d]=keqWorkArray[%5d]*Kwork[%5d];\n",fromKeqStepList[j].stepIdx,j,
                                                                    fromKeqStepList[j].fwdStepIdx);
    }
  fprintf(fptr,"}\n");
  fprintf(fptr,"\n");
  fprintf(fptr,"\n");
}


void rate_const::updateThirdBodyRxn(const double C[])
{
  int j,k;
  double Cmult;
  for(j=0; j<nThirdBodyRxn; j++)
    {
      //printf("3rd body reaction %d of %d:\n",j,nThirdBodyRxn);
      //printf("   nEnhanced    = %d\n", thirdBodyRxnList[j].nEnhanced);
      //printf("   fwd step id  = %d\n", thirdBodyRxnList[j].fwdStepIdx);
      //printf("   rev step id  = %d\n", thirdBodyRxnList[j].revStepIdx);

      Cmult=Csum;
      for(k=0; k<thirdBodyRxnList[j].nEnhanced; k++)
	{
	  //printf("      %d enhanced spc id = %d, eff = %g\n",
	  //	 k,thirdBodyRxnList[j].etbSpcIdx[k],
	  //	 thirdBodyRxnList[j].etbSpcEff[k]);
	  //fflush(stdout);

	  Cmult+=C[thirdBodyRxnList[j].etbSpcIdx[k]]*
	    thirdBodyRxnList[j].etbSpcEff[k];
	}
      Kwork[thirdBodyRxnList[j].fwdStepIdx]*=Cmult;
      if(likely(thirdBodyRxnList[j].revStepIdx >= 0))
	{Kwork[thirdBodyRxnList[j].revStepIdx]*=Cmult;}
    }
}

void rate_const::updateFalloffRxn(const double C[])
{
  int j,k;
  double Cmult,Klow,Pr,log_10_Pr,Pcorr,Fcenter,fTerm,nTerm;
  if(Tchanged) {
    for(j=0; j<nFalloffRxn; j++) {
      falloffWorkArray[j] = falloffRxnList[j].param[0]+
                            falloffRxnList[j].param[1]*log_e_Tcurrent-
                            falloffRxnList[j].param[2]*invTcurrent;
    }
    fast_vec_exp(falloffWorkArray, nFalloffRxn);
  }
  for(j=0; j<nFalloffRxn; j++) {

    if(falloffRxnList[j].falloffSpcIdx >= 0) {
      // single falloff species
      Cmult=C[falloffRxnList[j].falloffSpcIdx];
    }
    else {
      // third-body species falloffRxnList[j].falloffSpcIdx == MIN_INT32
      Cmult=Csum;
      for(k=0; k<falloffRxnList[j].nEnhanced; k++) {
	Cmult+=C[falloffRxnList[j].etbSpcIdx[k]]*
               falloffRxnList[j].etbSpcEff[k];
      }
    }

    Klow = falloffWorkArray[j];
    Pr = Klow*Cmult/Kwork[falloffRxnList[j].fwdStepIdx];
    if(Pr < 1.0e-300) {
      Pr = 1.0e-300; // ck SMALL constant
    }
    log_10_Pr=log10(Pr);

    fTerm=1.0; // default is Lindemann

    if(falloffRxnList[j].falloffType == TROE_THREE_PARAMS ||
       falloffRxnList[j].falloffType == TROE_FOUR_PARAMS) {

      // Troe 3 and 4-parameter fits
      Fcenter = 0.0;
      if(falloffRxnList[j].param[4]!=0) {
        Fcenter += (1.0-falloffRxnList[j].param[3])
                   *exp(-Tcurrent/falloffRxnList[j].param[4]);
      }
      if(falloffRxnList[j].param[5]!=0) {
        Fcenter+=falloffRxnList[j].param[3]
                 *exp(-Tcurrent/falloffRxnList[j].param[5]);
      }

      // Below are the special TROE alterations that were present
      // in JY Chen's version of chemkin II.  They are no longer used
      // because in one case, when alpha is less than zero, is actually used
      // in the full TROE form for the reaction C2H4+H(+M)<=>C2H5(+M)
      // reported by Miller and Klippenstein, Phys Chem Chem Phys, vol 6,
      // 1192-1202, 2004.
      //
      //if(falloffRxnList[j].param[4] < 0.0) {
      //  // Fcenter = T***
      //  Fcenter = -falloffRxnList[j].param[4];
      //}
      //
      //if(falloffRxnList[j].param[3] < 0.0) {
      //  // Fcenter = |alpha| + T*(T***)
      //	Fcenter = fabs(falloffRxnList[j].param[3])
      //            +falloffRxnList[j].param[4]*Tcurrent;
      //}

      if(falloffRxnList[j].falloffType ==  TROE_FOUR_PARAMS) {
        // 4-parameter Troe
        Fcenter += exp(-falloffRxnList[j].param[6]*invTcurrent);
      }

      // use original formulation
      if(Fcenter < 1.0e-300) {
        Fcenter = 1.0e-300;
      }
      fTerm=log10(Fcenter);
      nTerm=0.75-1.27*fTerm;
      log_10_Pr-=(0.4+0.67*fTerm);                // log10(Pr) + c
      log_10_Pr=log_10_Pr/(nTerm-0.14*log_10_Pr); // d = 0.14
      log_10_Pr*=log_10_Pr;
      fTerm/=(1.0+log_10_Pr);
      fTerm=pow(10.0,fTerm);
    // end if Troe 3 and 4 parameter falloff reactions
    } else if(falloffRxnList[j].falloffType == SRI) {

      // Note falloffRxnList[j].param[0-2] are the Klow/Khigh arrhenius
      // parameters.
      //
      const double a = falloffRxnList[j].param[3];
      const double b = falloffRxnList[j].param[4];
      const double inv_c = falloffRxnList[j].param[5];
      const double x_power = 1.0/(1.0+log_10_Pr*log_10_Pr);

      // Standard 3-term SRI definition
      //   F = (a*exp(-b/T) + exp(-T/c))**X
      fTerm = a*exp(-b*invTcurrent);
      if(inv_c > 0) {
        fTerm += exp(-Tcurrent*inv_c);
      }
      fTerm = pow(fTerm, x_power);

      // Auxillary 4 and 5-term SRI definitions
      //   F = d*(a*exp(-b/T) + exp(-T/c))**X         (4-term)
      //   F = d*(a*exp(-b/T) + exp(-T/c))**X * T**e  (5-term)
      // Note that the 4-term SRI function is not supported by Cantera
      // or Chemkin II.
      if(falloffRxnList[j].param.size() >= 7) {
        fTerm *= falloffRxnList[j].param[6];  // pre-multiplier 'd'
      }
      if(falloffRxnList[j].param.size() == 8) {
        fTerm *= pow(Tcurrent, falloffRxnList[j].param[7]); // multiplier T**e
      }
    }

    Pcorr = fTerm*Pr/(1.0+Pr);

    Kwork[falloffRxnList[j].fwdStepIdx]*=Pcorr;
    if(likely(falloffRxnList[j].revStepIdx >= 0))
      {Kwork[falloffRxnList[j].revStepIdx]*=Pcorr;}
  }
}

void rate_const::getKrxn(info_net &netobj, double Kfwd[],
			     double Krev[])
{
  int j;

  for(j=0; j<netobj.getNumRxn(); j++)
    {Kfwd[j]=Krev[j]=0.0;}
  for(j=0; j<nStep; j++)
    {
      if(netobj.getRxnDirOfStep(j)==1)
	{Kfwd[netobj.getRxnIdxOfStep(j)]=Kwork[j];}
      else
	{Krev[netobj.getRxnIdxOfStep(j)]=Kwork[j];}
    }
}

int rate_const::getThirdBodyEff(ckr::CKReader &ckrobj, int rxnId,
				    vector <int> &spcId,
				    vector <double> &spcEff)
{
  int id,j,k;
  int nEnh=ckrobj.reactions[rxnId].e3b.size();
  string spcName;
  map<string, double>::iterator iter;

  k=0;
  iter=ckrobj.reactions[rxnId].e3b.begin();
  for(j=0; j<nEnh; j++)
    {
      spcName=(*iter).first;
      id=spcIdxOfString(ckrobj,spcName);
      if(id < 0)
	{
	  printf("WARNING: reaction %d enhanced third body species %s\n",
		 rxnId+1,spcName.c_str());
	  printf("	   not found in mechanism - ignoring it\n");
	}
      else
	{
	  spcId[k]=id;
	  spcEff[k]=(*iter).second-1.0;
	  ++k;
	}
      (*iter++);
    }
  //spcId.resize(k);
  //spcEff.resize(k);
  return k;
}


int rate_const::spcIdxOfString(ckr::CKReader &ckrobj, string spcName)
{
  int j;

  for(j=0; j<nSpc; j++)
    {
      if(ckrobj.species[j].name == spcName)
	{return j;}
    }
  return MIN_INT32;
}
void rate_const::print()
{
  printf("# From void rate_const::print()\n");
  printf("#\n");
  printf("# number of chemical reaction steps (1-way)    : %d\n",nStep);
  printf("# number of steps with Arrhenius T-func        : %d\n",nArrheniusStep);
  printf("# number of distinct Arrhenius T-func          : %d\n",
	 nDistinctArrhenius);
  printf("# number of steps with Landau-Teller T-func    : %d\n",nLandauTellerStep);
  printf("# number of steps with Krev from Keq T-func    : %d\n",nFromKeqStep);
  printf("# number of reactions with third-body          : %d\n",
	 nThirdBodyRxn);
  printf("# number of reactions with falloff             : %d\n",
	 nFalloffRxn);


}
int isSameArrheniusTol(arrheniusSortElem x, arrheniusSortElem y)
{
  if(fabs(x.A - y.A) > ARRHENIUS_RTOL*fabs(0.5*(x.A+y.A)) &&
     fabs(x.A - y.A) > ARRHENIUS_ATOL)
    {return 0;}
  if(fabs(x.Tpow - y.Tpow) > ARRHENIUS_RTOL*fabs(0.5*(x.Tpow+y.Tpow)) &&
     fabs(x.Tpow - y.Tpow) > ARRHENIUS_ATOL)
    {return 0;}
  if(fabs(x.Tact - y.Tact) > ARRHENIUS_RTOL*fabs(0.5*(x.Tact+y.Tact)) &&
      fabs(x.Tact - y.Tact) > ARRHENIUS_ATOL)
    {return 0;}

  return 1;
}
int compareArrhenius(const void *x, const void *y)
{
  arrheniusSortElem *xptr = (arrheniusSortElem *)x;
  arrheniusSortElem *yptr = (arrheniusSortElem *)y;

  if(xptr->A < yptr->A)
    {return -1;}
  else if(xptr->A > yptr->A)
    {return 1;}

  // note a larger activation temperature yields a smaller rate constant
  if(xptr->Tact > yptr->Tact)
    {return -1;}
  else if(xptr->Tact < yptr->Tact)
    {return 1;}

  if(xptr->Tpow < yptr->Tpow)
    {return -1;}
  else if(xptr->Tpow > yptr->Tpow)
    {return 1;}

  return 0;
 }

int compareArrheniusT1000(const void *x, const void *y)
{
  arrheniusSortElem *xptr = (arrheniusSortElem *)x;
  arrheniusSortElem *yptr = (arrheniusSortElem *)y;

  double Kx=exp(log(xptr->A)+xptr->Tpow*log(1000.0)-xptr->Tact/1000.0);
  double Ky=exp(log(yptr->A)+yptr->Tpow*log(1000.0)-yptr->Tact/1000.0);

  if(Kx < Ky)
    {return -1;}
  else if(Kx > Ky)
    {return 1;}
  return 0;
}

void rate_const::setPLogInterpolationStepList(ckr::CKReader &ckrobj,
                                  info_net &netobj)
{
  const int nRxn=netobj.getNumRxn();
  int step_id,num_pts;
  std::vector<double> pres;
  std::vector<double> Afactor;
  std::vector<double> Tpow;
  std::vector<double> Tact;

  for(int j=0; j<nRxn; ++j) {
    if(ckrobj.reactions[j].kf.type == ckr::PLogInterpolation) {
      step_id = netobj.getStepIdxOfRxn(j,1);

      num_pts = ckrobj.reactions[j].kf.pres_pts_plog.size();
      pres    = ckrobj.reactions[j].kf.pres_pts_plog;
      Afactor = ckrobj.reactions[j].kf.A_plog;
      Tpow    = ckrobj.reactions[j].kf.n_plog;
      Tact    = ckrobj.reactions[j].kf.E_plog;

      for(int k=0; k<num_pts; ++k) {
        pres[k]*=P_ATM;    // convert mech file PLOG pressure [atm] -> [Pa]
        Tact[k]*=convertE; // convert mech file activation energy to [K]
        Afactor[k]*=pow(convertC,
                        (1.0-netobj.getRealOrderOfStep(step_id)));
        // there is no 3rd body units correction since PLOG and 3rd body
        // reactions are treated as being mutually exclusive
      }
      PLogReaction current_reaction(PLogReaction(j,
                                                 step_id,
                                                 num_pts,
                                                 pres.data(),
                                                 Afactor.data(),
                                                 Tpow.data(),
                                                 Tact.data(),
                                                 false)); // use_extrapolate

      plogInterpolationStepList.push_back(current_reaction);
    }
  }
}

void rate_const::updatePLogInterpolationStep(const double pressure,
                                             const double log_e_pressure)
{
  for(int j=0; j<nPLogInterpolationStep; ++j) {

   Kwork[plogInterpolationStepList[j].step_index()] =
     plogInterpolationStepList[j].GetRateCoefficientFromTP(Tcurrent,
                                                           invTcurrent,
                                                           log_e_Tcurrent,
                                                           pressure,
                                                           log_e_pressure);

  }
}
} // namespace zerork
