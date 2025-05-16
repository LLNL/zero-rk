#ifndef ZERORK_MECHANISM_H
#define ZERORK_MECHANISM_H

#include <string>
#include "../CKconverter/CKReader.h"
#include "element.h"
#include "species.h"
#include "info_net.h"
#include "perf_net.h"
#include "rate_const.h"
#include "constants.h"
#include "nasa_poly.h"
#include "atomicMassDB.h"
#include "external_funcs.h"
#include "non_integer_reaction_network.h"
#include "utilities.h"

namespace zerork {

class mechanism
{
 public:
  mechanism(const char *mechFileName,
		const char *thermFileName,
		const char *convertFileName,
    int verbosity_inp = 1);
  virtual ~mechanism();

  int getIdxFromName(const char *nm);
  int getNumElements() const {return nElm;}
  std::map<std::string, double> getElementInfo() const;
  std::map<std::string, std::map<std::string, int> > getSpeciesElementInfo() const;
  int getNumSpecies() const {return nSpc;}
  int getNumReactions() const {return nRxn;}
  int getNumSteps() const {return nStep;}
  int getRxnIdxOfStep(const int stepId) const
  {return infoNet->getRxnIdxOfStep(stepId);}
  void getStepIdxOfRxn(const int rxnId, int *fwdId, int *revId)
  {
    (*fwdId) = infoNet->getStepIdxOfRxn(rxnId,1);
    (*revId) = infoNet->getStepIdxOfRxn(rxnId,-1);
  }
  int getStepIdxOfRxn(int id, int dir) const //alternate version
        {return infoNet->getStepIdxOfRxn(id,dir);}
  int getOrderOfStep(int id) const {return infoNet->getOrderOfStep(id);}
  int getNumProductsOfStep(int id) const {return infoNet->getNumProductsOfStep(id);}
  int getSpecIdxOfStepReactant(const int stepId, const int reacId) const
        {return infoNet->getSpecIdxOfStepReactant(stepId,reacId);}
  int getSpecIdxOfStepProduct(const int stepId, const int prodId) const
        {return infoNet->getSpecIdxOfStepProduct(stepId,prodId);}
  int getTotalReactants() const {return infoNet->getTotalReactants();}
  int getTotalProducts() const {return infoNet->getTotalProducts();}
  int getMaxProductsInStep() const {return infoNet->getMaxProductsInStep();}
  int getMaxReactantsInStep() const {return infoNet->getMaxReactantsInStep();}

  const char * getSpeciesName(const int idx) const;
  const char * getReactionName(const int idx) const
  {return rxnDefinition[idx].c_str();}
  void getReactionNameDirOfStep(const int idx, std::string *str) const;
  string getReactionString(const int idx) const
  {return rxnDefinition[idx];}

  double getGasConstant() const {return Ru;}
  double getAvogadroNumber() const {return KAvogadroNumber;}

  double getMolWtMixFromX(const double x[]) const;
  double getMolWtMixFromY(const double y[]) const;
  void getMolWtSpc(double mwCopy[]) const;

  double getDensityFromTPY(const double T, const double P,
			   const double y[]) const;
  double getDensityFromTPX(const double T, const double P,
			   const double x[]) const;
  double getPressureFromTVY(const double T, const double v,
			    const double y[]) const;
  //double getPressureFromTC(const double T, const double conc[]) const;

  void getYfromX(const double x[], double y[]) const;
  void getXfromY(const double y[], double x[]) const;
  void getXfromC(const double conc[], double x[]) const;
  void getCfromVY(const double v, const double y[], double c[]) const;

  void getThermoCoeffs(double coeffs[]) const;

  double getMassEnthalpyFromTY(const double T, const double y[],
			       double hSpc[]) const;
  double getMassEnthalpyFromTY(const double T, const double y[]) const;
  double getMassIntEnergyFromTY(const double T, const double y[],
				double uSpc[]) const;
  double getMassIntEnergyFromTY(const double T, const double y[]) const;
  void getEnthalpy_RT(const double T, double h_RT[]) const;
  void getIntEnergy_RT(const double T, double u_RT[]) const;
  void getCp_R_Enthalpy_RT(const double T,
                           double cp_R[],
                           double h_RT[]) const;
  void getCv_R_IntEnergy_RT(const double T,
                            double cv_R[],
                            double u_RT[]) const;
  double getMassCpFromTY(const double T, const double y[],
			 double cpSpc[]) const;
  double getMassCpFromTY(const double T, const double y[]) const;
  double getMassCvFromTY(const double T, const double y[],
			 double cvSpc[]) const;
  double getMassCvFromTY(const double T, const double y[]) const;
  double getMolarCvFromTC(const double T, const double c[],
			 double cvSpc[]) const;
  double getMolarCvFromTC(const double T, const double c[]) const;
  void getNonDimGibbsFromT(const double T, double G_RT[]) const;
  double getTemperatureFromEY(const double E, const double y[], const double temp_guess) const;
  double getTemperatureFromHY(const double H, const double y[], const double temp_guess) const;

  void getKrxnFromTC(const double T, const double C[], double Kfwd[],
		     double Krev[]);
  void getReactionRates(const double T, const double C[],
			double netOut[], double createOut[],
			double destroyOut[], double stepOut[]);
  void getReactionRatesFromTCM(const double T,
                               const double C[],
                               const double C_mix,
			       double netOut[],
                               double createOut[],
			       double destroyOut[],
                               double stepOut[]);
  void getNetReactionRates(const double T,
                           const double C[],
                           double netOut[]);
  void getReactionRates_perturbROP(const double T, const double C[],
		                   const double perturbMult[],
			           double netOut[], double createOut[],
			           double destroyOut[], double stepOut[]);

  // Compute the reaction rates and the rate-of-progress of each step after
  // applying the step limiter to rate coefficient.  The step limiter is
  // applied as follows (K is the rate coefficient, K_lim is the limited rate
  // coefficient used to compute the rate of progress of each step, and L is
  // the limit value stored in the step_limiter array):
  //
  //    K_lim = K * [L/(K + L)]  note K_lim <= L
  //
  // the relative error using the limiter is (K/L)/((K/L) + 1) so when
  // K is 1/100 of L the relative error is 0.01/1.01 or ~1%
  //
  void getReactionRatesLimiter(const double T,
                               const double C[],
		               const double step_limiter[],
			       double netOut[],
                               double createOut[],
			       double destroyOut[],
                               double stepOut[]);

  void getReactionRatesLimiter_perturbROP(const double T,
                                          const double C[],
                                          const double step_limiter[],
                                          const double perturbMult[],
                                          double netOut[],
                                          double createOut[],
                                          double destroyOut[],
                                          double stepOut[]);

  // these may need to be inherited
  double getMolarAtomicOxygenRemainder(const double x[]) const;
  void getMolarIdealLeanExhaust(const double xInit[],
			    double xFinal[]) const;
  void getMolarIdealExhaust(const double xInit[],
			    double xFinal[]) const;

  double getProgressEquivalenceRatio(const double y[]) const;
  void getModifiedEquivalenceRatios(const double y[],
                                    double &smallmxf,
                                    double &bigmxf) const;
  double getCHValue(const double y[]) const;
  void getOxygenAtomCount(int numO[]) const;
  void getCarbonAtomCount(int numC[]) const;
  void getHydrogenAtomCount(int numH[]) const;
  // TODO: resolve duplicate functionality in the atomic count functions
  //       currently leaving both sets to facilitate merge
  void getSpeciesHydrogenCount(int num_atoms[]) const;
  void getSpeciesNitrogenCount(int num_atoms[]) const;
  void getSpeciesCarbonCount(int num_atoms[]) const;
  void getSpeciesOxygenCount(int num_atoms[]) const;
  void getSpeciesArgonCount(int num_atoms[]) const;
  void getSpeciesHeliumCount(int num_atoms[]) const;

  // load external func lib if called for
  void initExternalFuncs();

  // multi-reactor functions
  void getMassCpFromTY_mr(const int nReactors, const double T[],
                          const double y[], double cpSpc[],
                          double cpReactors[]) const;
  void getMassCvFromTY_mr(const int nReactors, const double T[],
                          const double y[], double cvSpc[],
                          double cvReactors[]) const;
  void getEnthalpy_RT_mr(const int nReactors, const double T[], double h_RT[]) const;
  void getIntEnergy_RT_mr(const int nReactors, const double T[],
                          double u_RT[]) const;
  void getCfromVY_mr(const int nReactors, const double v[],
                     const double y[], double c[]) const;
  void getDensityFromTPY_mr(const int nReactors, const double *T,
                            const double *P, const double y[],
                            double *dens) const;
  void getMolWtMixFromY_mr(const int nReactors, const double y[],
                           double *mwMix) const;

  int isThirdBodyReaction(const int rxn_id) const
  {return ((infoNet->getThirdBodyFlagOfReaction(rxn_id) == 0) ? 0 : 1);}

  int isFalloffReaction(const int rxn_id) const
  {return ((infoNet->getFalloffFlagOfReaction(rxn_id) == 0) ? 0 : 1);}

  int isThirdBodyFalloffReaction(const int rxn_id) const
  {return ((infoNet->getFalloffFlagOfReaction(rxn_id) == -1) ? 1 : 0);}

  int isReversibleReaction(const int rxn_id) const
  {return ((infoNet->getReversibleFlagOfReaction(rxn_id) == 0) ? 0 : 1);}

  int getNumEnhancedSpeciesOfStep(const int id) const
  {return infoNet->getNumEnhancedSpeciesOfStep(id);}

  int getEnhancementFactorsOfStep(const int id,
                                  std::vector<int> *species_id,
                                  std::vector<double> *alpha)
  {return infoNet->getEnhancementFactorsOfStep(id,species_id,alpha);}

  void generateCircosFilesTPY(const double T, const double P,
                              const double y[],
                              const int num_reacs,
                              const int exclude_unimolecular,
                              const char *filename_prefix
                              );
  NonIntegerReactionNetwork * getNonIntegerReactionNetwork()
  {return &non_integer_network_;}

 protected:
  void buildReactionString(const int idx,
                           string &str);

  int nSpc;    // # of species
  // performance arrays
  double *molWt;
  double *invMolWt;
  double *RuInvMolWt;

  void initialize_ptrs(ckr::CKReader *ckrobj);//allow subclass to use different ptrs
  nasa_poly_group *thermo;
  rate_const *Kconst;
  perf_net *perfNet;
  info_net *infoNet;
  NonIntegerReactionNetwork non_integer_network_;
  string *rxnDefinition;

  //File name for mechanism, thermo, parser output
  std::string mechFileStr;
  std::string thermFileStr;
  std::string convertFileStr;
  int verbosity;

 private:
  void build_mechanism(ckr::CKReader *ckrobj);

  // sizes and constants
  int nElm;    // # of elements
  int nRxn;    // # of 2-way reactions
  int nStep;
  double Ru;   // [J/kmol-K] universal gas constant

  // informational, organizational data
  species *speciesList;
  std::map<std::string, double> elementInfo;
  std::map<std::string, std::map<std::string, int> > speciesElementInfo;

  // reaction rate workspace
  double *createRateWorkspace;
  double *destroyRateWorkspace;
  double *stepROPWorkspace;

  // handles for loading funcs from external lib
  void * externalFuncLibHandle;
  external_func_check_t ex_func_check; // test out loaded library,
                                       // args are nsp nreacsteps
  external_func_rates_t ex_func_calc_rates;
  external_func_arrh_t ex_func_calc_arrh;
  external_func_keq_t ex_func_calc_keq;
};

} // namespace zerork

#endif
