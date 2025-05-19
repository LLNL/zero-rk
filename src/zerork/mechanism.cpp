#include <stdlib.h>
#include <string.h> // for memcpy
#include <assert.h>
#ifndef _WIN32
#include <dlfcn.h> // for loading external func lib
#endif

#include <string>
#include <algorithm> //for std::sort
#include <exception>

#include "mechanism.h"

namespace zerork {


mechanism::mechanism(const char *mechFileName,
                     const char *thermFileName,
                     const char *convertFileName,
                     int verbosity_inp)
    :
  mechFileStr(mechFileName),
  thermFileStr(thermFileName),
  convertFileStr(convertFileName),
  verbosity(verbosity_inp)
{
  //Note: if strlen(convertFileName) == 0 we will not print parser output
  //      unless an error happens.  Then we will print parser output
  //      only if verbosity > 0.
  //      This is our strategy for multiple parsers on MPI programs.
  ckr::CKReader ckrobj;
  if(!(ckrobj.read(mechFileStr,thermFileStr,convertFileStr)))
  {
    if(convertFileStr.length() != 0)
    {
      printf("ERROR: could not parse mech (%s) and thermo (%s) files,\n",
             mechFileStr.c_str(),thermFileStr.c_str());
      printf("       check converter log file %s\n",convertFileStr.c_str());
    }
    else if(verbosity > 0)
    {
      //Retry parse with non-null output
      convertFileStr = mechFileStr + ".out";
      ckr::CKReader ckrobj2;
      if(!(ckrobj2.read(mechFileStr,thermFileStr,convertFileStr)))
      {
        printf("ERROR: could not parse mech (%s) and thermo (%s) files,\n",
               mechFileStr.c_str(),thermFileStr.c_str());
        printf("       check converter log file %s\n",convertFileStr.c_str());
      }
    }
    fflush(stdout);
#ifdef ZERORK_MECHANISM_EXIT_FAIL
    exit(1);
#else
    throw std::runtime_error("zerork::mechanism::mechanism() failed to parse mechanism.");
#endif
  }

  //This is our flag to not output parser warnings from rate_const
  if(verbosity == 0 || convertFileStr.length() == 0)
  {
    ckrobj.verbose = false;
  }
  build_mechanism(&ckrobj);
}

mechanism::~mechanism()
{
  delete perfNet;
  delete Kconst;
  delete infoNet;
  delete [] rxnDefinition;
  delete [] speciesList;
  delete [] molWt;
  delete [] invMolWt;
  delete [] RuInvMolWt;

  delete [] createRateWorkspace;
  delete [] destroyRateWorkspace;
  delete [] stepROPWorkspace;

  delete thermo;
  //delete rxnNet;
#ifndef _WIN32
  if(externalFuncLibHandle) dlclose(externalFuncLibHandle);
#endif
}

// build_mechanism(...)
//
// This function contains all the constructor steps. This allows the
// bulk of the constructor steps to be located in one location.
void mechanism::build_mechanism(ckr::CKReader *ckrobj)
{
  int nConst;
  int j,k;
  vector <int> constCount;
  vector <element> constList;

  assert(ckrobj != NULL);

  // set the gas constant
  Ru=NIST_RU;

  // create and fill the element list
  nElm=static_cast<int>(ckrobj->elements.size());
  element dummyElement;
  for(j=0; j<nElm; ++j) {
    bool flag = dummyElement.setElement(ckrobj->elements[j].name.c_str());
    elementInfo[dummyElement.getSymbol_c_str()] = dummyElement.getMass();
    assert(("Unable to set element", flag));
  }

  // create and fill the species sized lists
  nSpc=static_cast<int>(ckrobj->species.size());
  nRxn=static_cast<int>(ckrobj->reactions.size());
  speciesList = new species[nSpc];
  molWt=new double[nSpc];
  invMolWt=new double[nSpc];
  RuInvMolWt=new double[nSpc];

  for(j=0; j<nSpc; j++)
  {
      // set Species.index to its position in the species list and the species
      // name map
      ckrobj->species[j].index=j;
      ckrobj->speciesData[ckrobj->species[j].name].index=j;

      nConst=ckrobj->species[j].elements.size(); // number of constituent elements
      constCount.resize(nConst);
      constList.resize(nConst);
      std::map<std::string, int> elemCountMap;
      for(k=0; k<nConst; k++)
      {
        constCount[k]=ckrobj->species[j].elements[k].number;
        bool flag = constList[k].setElement(ckrobj->species[j].elements[k].name.c_str());
        elemCountMap[constList[k].getSymbol_c_str()] = constCount[k];
        assert(("Unable to set element",flag));
      }
      speciesList[j].setSpecies(j,ckrobj->species[j].name.c_str(),constCount,
				constList);
      speciesElementInfo[speciesList[j].getName_c_str()] = elemCountMap;

      molWt[j]=speciesList[j].getMolecularWeight();
      invMolWt[j]=1.0/molWt[j];
      RuInvMolWt[j]=Ru/molWt[j];
  }



  initialize_ptrs(ckrobj);
  non_integer_network_ = infoNet->getNonIntegerReactionNetwork();
  //printf("# INFO: Number of non-integer reactions = %d\n",
  //       non_integer_network.GetNumNonIntegerReactions());
  //printf("# INFO: Number of non-integer steps     = %d\n",
  //       non_integer_network.GetNumNonIntegerSteps());


  nStep=infoNet->getNumSteps();
  // workspace for calculating the reaction rates
  createRateWorkspace  = new double[nSpc];
  destroyRateWorkspace = new double[nSpc];
  stepROPWorkspace     = new double[nStep];

  rxnDefinition = new string[nRxn];
  for(j=0; j<nRxn; j++) {
    if(non_integer_network_.HasReaction(j)) {
      rxnDefinition[j] = non_integer_network_.GetNameOfReaction(j);
    } else {
      buildReactionString(j,rxnDefinition[j]);
    }
    //cout << "rxnDefinition[" << j<< "] " << rxnDefinition[j] << endl;
  }

  initExternalFuncs();

}

void mechanism::initialize_ptrs(ckr::CKReader *ckrobj)
{
  int j,k;
  int nLoCoef,nHiCoef;
  double *Tlo, *Thi, *coef;

  // fill and create the data for the thermodynamics object
  coef = new double[nSpc*NUM_THERMO_POLY_D5R2];
  Tlo  = new double[nSpc];
  Thi  = new double[nSpc];

  for(j=0; j<nSpc; j++)
  {
    Tlo[j]=ckrobj->species[j].tlow;
    Thi[j]=ckrobj->species[j].thigh;
    coef[j*NUM_THERMO_POLY_D5R2+0]=ckrobj->species[j].tmid;

    nLoCoef=ckrobj->species[j].lowCoeffs.size();
    assert(nLoCoef == NUM_COEF_RANGE);
    for(k=0; k<NUM_COEF_RANGE; k++)
    {
      coef[j*NUM_THERMO_POLY_D5R2+k+1]= ckrobj->species[j].lowCoeffs[k];
    }
    nHiCoef=ckrobj->species[j].highCoeffs.size();
    assert(nHiCoef == NUM_COEF_RANGE);
    for(k=0; k<NUM_COEF_RANGE; k++)
    {
      coef[j*NUM_THERMO_POLY_D5R2+k+1+NUM_COEF_RANGE]=
          ckrobj->species[j].highCoeffs[k];
    }
  }

  thermo = new nasa_poly_group(nSpc,coef,Tlo,Thi);
  infoNet = new info_net(ckrobj);
  Kconst = new rate_const(ckrobj,infoNet,thermo);
  perfNet = new perf_net(*infoNet,*Kconst);

  delete [] Thi;
  delete [] Tlo;
  delete [] coef;
}

int mechanism::getIdxFromName(const char *nm)
{
  int idx;
  for(idx=0; idx<nSpc; idx++) {
    if(strcmp(nm,speciesList[idx].getName_c_str())==0) {
      return idx;
    }
  }
  return -1;
}


std::map<std::string, double> mechanism::getElementInfo() const {
   return elementInfo;
}

std::map<std::string, std::map<std::string, int> > mechanism::getSpeciesElementInfo() const {
   return speciesElementInfo;
}

const char * mechanism::getSpeciesName(const int idx) const
{return speciesList[idx].getName_c_str();}

void mechanism::getReactionNameDirOfStep(const int idx,
                                         std::string *str) const
{
  int reaction_id = getRxnIdxOfStep(idx);
  *str = std::string(getReactionName(reaction_id));
  if(getStepIdxOfRxn(reaction_id,1) == idx) {
    (*str)+=" (fwd)";
  } else {
    (*str)+=" (rev)";
  }
}


double mechanism::getMolWtMixFromX(const double x[]) const
{
  double mwSum=0.0;

  for(int j=0; j<nSpc; j++)
    {mwSum+=x[j]*molWt[j];}
  return mwSum;
}

double mechanism::getMolWtMixFromY(const double y[]) const
{
  double mwSum=0.0;

  for(int j=0; j<nSpc; j++)
    {mwSum+=y[j]*invMolWt[j];}
  return 1.0/mwSum;
}

void mechanism::getMolWtSpc(double mwCopy[]) const
{memcpy(mwCopy,molWt,sizeof(double)*nSpc);}

double mechanism::getDensityFromTPY(const double T, const double P,
					const double y[]) const
{
  double molWtMix=getMolWtMixFromY(y);
  return P*molWtMix/(Ru*T);
}

double mechanism::getDensityFromTPX(const double T, const double P,
					const double x[]) const
{
  double molWtMix=getMolWtMixFromX(x);
  return P*molWtMix/(Ru*T);
}
double mechanism::getPressureFromTVY(const double T, const double v,
					 const double y[]) const
{
  double molWtMix=getMolWtMixFromY(y);
  return Ru*T/(v*molWtMix);
}

void mechanism::getYfromX(const double x[], double y[]) const
{
  double invMolWtMix=1.0/getMolWtMixFromX(x);
  for(int j=0; j<nSpc; j++)
    {y[j]=invMolWtMix*x[j]*molWt[j];}
}
void mechanism::getXfromY(const double y[], double x[]) const
{
  double molWtMix=getMolWtMixFromY(y);
  for(int j=0; j<nSpc; j++)
    {x[j]=y[j]*molWtMix*invMolWt[j];}
}

void mechanism::getXfromC(const double conc[], double x[]) const
{
  double Csum=0.0;
  for(int j=0; j<nSpc; j++)
    {Csum+=conc[j];}
  Csum=1.0/Csum;
  for(int j=0; j<nSpc; j++)
    {x[j]=conc[j]*Csum;}
}


void mechanism::getCfromVY(const double v, const double y[],
			       double c[]) const
{
  double dens=1.0/v;
  for(int j=0; j<nSpc; j++)
    {c[j]=dens*invMolWt[j]*y[j];}
}

void mechanism::getThermoCoeffs(double coeffs[]) const
{
  thermo->getThermoCoeffs(coeffs);
}

double mechanism::getMassEnthalpyFromTY(const double T, const double y[],
					    double hSpc[]) const
{
  double hMix=0.0;
  thermo->getH_RT(T,hSpc);
  for(int j=0; j<nSpc; j++)
    {
      hSpc[j]*=RuInvMolWt[j]*T;
      hMix+=hSpc[j]*y[j];
    }
  return hMix;
}

double mechanism::getMassEnthalpyFromTY(const double T, const double y[]) const
{
  std::vector<double> hSpc(nSpc);
  return getMassEnthalpyFromTY(T, y, &hSpc[0]);
}

double mechanism::getMassIntEnergyFromTY(const double T, const double y[],
					    double uSpc[]) const
{
  double uMix=0.0;
  thermo->getH_RT(T,uSpc);
  for(int j=0; j<nSpc; j++)
    {
      uSpc[j]-=1.0;
      uSpc[j]*=RuInvMolWt[j]*T;
      uMix+=uSpc[j]*y[j];
    }
  return uMix;
}

double mechanism::getMassIntEnergyFromTY(const double T, const double y[]) const
{
  std::vector<double> eSpc(nSpc);
  return getMassIntEnergyFromTY(T, y, &eSpc[0]);
}


void mechanism::getEnthalpy_RT(const double T, double h_RT[]) const
{
  thermo->getH_RT(T,h_RT);
}

void mechanism::getIntEnergy_RT(const double T, double u_RT[]) const
{
  thermo->getH_RT(T,u_RT);
  for(int j=0; j<nSpc; j++)
    {u_RT[j]-=1.0;}
}

void  mechanism::getCp_R_Enthalpy_RT(const double T,
                                         double cp_R[],
                                         double h_RT[]) const
{
  thermo->getCp_R(T,&cp_R[0]);
  thermo->getH_RT(T,&h_RT[0]);
}
void  mechanism::getCv_R_IntEnergy_RT(const double T,
                                          double cv_R[],
                                          double u_RT[]) const
{
  thermo->getCp_R(T,&cv_R[0]);
  thermo->getH_RT(T,&u_RT[0]);
  for(int j=0; j<nSpc; ) {
    cv_R[j]-=1.0;
    u_RT[j]-=1.0;
    ++j;
  }
}

double mechanism::getMassCpFromTY(const double T, const double y[],
				      double cpSpc[]) const
{
  double cpMix=0.0;
  thermo->getCp_R(T,cpSpc);
  for(int j=0; j<nSpc; j++)
  {
    cpSpc[j]*=RuInvMolWt[j];
    cpMix+=cpSpc[j]*y[j];
  }

  return cpMix;
}

double mechanism::getMassCpFromTY(const double T, const double y[]) const
{
  std::vector<double> cpSpc(nSpc);
  return getMassCpFromTY(T,y,&cpSpc[0]);
}

double mechanism::getMassCvFromTY(const double T, const double y[],
				      double cvSpc[]) const
{
  double cvMix=0.0;
  thermo->getCp_R(T,&cvSpc[0]);
  for(int j=0; j<nSpc; j++)
    {
      cvSpc[j]-=1.0;  // Cv = Cp
      cvSpc[j]*=RuInvMolWt[j];
      cvMix+=cvSpc[j]*y[j];
    }

  return cvMix;
}

double mechanism::getMassCvFromTY(const double T, const double y[]) const
{
  std::vector<double> cvSpc(nSpc);
  return getMassCvFromTY(T,y,&cvSpc[0]);
}

double mechanism::getMolarCvFromTC(const double T, const double c[],
				       double cvSpc[]) const
{
  double cvMix=0.0;
  double csum=0.0;
  double Ru=getGasConstant();

  thermo->getCp_R(T,&cvSpc[0]);
  for(int j=0; j<nSpc; j++)
    {
      cvSpc[j]-=1.0;  // Cv = Cp
      cvSpc[j]*=Ru;
      cvMix+=cvSpc[j]*c[j];
      csum+=c[j];
    }

  return cvMix/csum;
}

double mechanism::getMolarCvFromTC(const double T, const double c[]) const
{
  std::vector<double> cvSpc(nSpc);
  return getMolarCvFromTC(T,c,&cvSpc[0]);
}

void mechanism::getNonDimGibbsFromT(const double T, double G_RT[]) const
{thermo->getG_RT(T,G_RT);}

double mechanism::getTemperatureFromEY(const double E, const double y[], const double temp_guess) const
{
    const int maximum_iterations = 400;
    const double tolerance  = 1.0e-6;
    const double min_temperature = 100;
    const double max_temperature = 5000;

    std::vector<double> tmp(nSpc);
    double min_E = getMassIntEnergyFromTY(min_temperature, &y[0], &tmp[0]);
    double max_E = getMassIntEnergyFromTY(max_temperature, &y[0], &tmp[0]);
    if (E < min_E) {
        //Extrapolate
        double cv = getMassCvFromTY(min_temperature, &y[0], &tmp[0]);
        return min_temperature - (min_E-E)/cv;
    }
    if (E > max_E) {
        //Extrapolate
        double cv = getMassCvFromTY(max_temperature, &y[0], &tmp[0]);
        return max_temperature - (max_E-E)/cv;
    }
    double temp_iter = temp_guess;
    if(temp_iter < min_temperature || temp_iter > max_temperature) {
      temp_iter = min_temperature + (max_temperature-min_temperature)/(max_E-min_E)*(E-min_E);
    }
    for(int i = 0; i < maximum_iterations; ++i) {
        double E_iter = getMassIntEnergyFromTY(temp_iter, &y[0], &tmp[0]);
        double cv = getMassCvFromTY(temp_iter, &y[0], &tmp[0]);
        double delta_temp = std::min(std::max( (E-E_iter)/cv, -100.), 100.);
        if(std::abs(delta_temp) < tolerance || temp_iter+delta_temp == temp_iter) {
           break;
        }
        temp_iter += delta_temp;
    }
    return temp_iter;
}

double mechanism::getTemperatureFromHY(const double H, const double y[], const double temp_guess) const
{
    const int maximum_iterations = 400;
    const double tolerance  = 1.0e-6;
    const double min_temperature = 100;
    const double max_temperature = 5000;

    std::vector<double> tmp(nSpc);
    double min_H = getMassEnthalpyFromTY(min_temperature, &y[0], &tmp[0]);
    double max_H = getMassEnthalpyFromTY(max_temperature, &y[0], &tmp[0]);
    if (H < min_H) {
        //Extrapolate
        double cp = getMassCpFromTY(min_temperature, &y[0], &tmp[0]);
        return min_temperature - (min_H-H)/cp;
    }
    if (H > max_H) {
        //Extrapolate
        double cp = getMassCvFromTY(max_temperature, &y[0], &tmp[0]);
        return max_temperature - (max_H-H)/cp;
    }
    double temp_iter = temp_guess;
    if(temp_iter < min_temperature || temp_iter > max_temperature) {
      temp_iter = min_temperature + (max_temperature-min_temperature)/(max_H-min_H)*(H-min_H);
    }
    for(int i = 0; i < maximum_iterations; ++i) {
        double H_iter = getMassEnthalpyFromTY(temp_iter, &y[0], &tmp[0]);
        double cp = getMassCvFromTY(temp_iter, &y[0], &tmp[0]);
        double delta_temp = std::min(std::max( (H-H_iter)/cp, -100.), 100.);
        if(std::abs(delta_temp) < tolerance || temp_iter+delta_temp == temp_iter) {
           break;
        }
        temp_iter += delta_temp;
    }
    return temp_iter;
}

void mechanism::getKrxnFromTC(const double T, const double C[],
				  double Kfwd[], double Krev[])
{
  Kconst->updateK(T,C);
  Kconst->getKrxn(*infoNet,&Kfwd[0],&Krev[0]);
}

void mechanism::getReactionRates(const double T, const double C[],
				     double netOut[], double createOut[],
				     double destroyOut[], double stepOut[])
{
  perfNet->calcRatesFromTC(T,&C[0],&netOut[0],&createOut[0],&destroyOut[0],
  			   &stepOut[0]);
  //perfNet->calcRatesFromTC_unroll16(T,&C[0],&netOut[0],&createOut[0],
  //				    &destroyOut[0],&stepOut[0]);
  //perfNet->calcRatesFromExplicit(T,&C[0],&netOut[0],&createOut[0],
  //				 &destroyOut[0],&stepOut[0]);
}

void mechanism::getReactionRatesFromTCM(const double T,
                                        const double C[],
                                        const double C_mix,
			                double netOut[],
                                        double createOut[],
			                double destroyOut[],
                                        double stepOut[])
{
  perfNet->calcRatesFromTCM(T,
                            &C[0],
                            C_mix,
                            &netOut[0],
                            &createOut[0],
                            &destroyOut[0],
  			    &stepOut[0]);
}

void mechanism::getNetReactionRates(const double T,
                                    const double C[],
                                    double netOut[])
{
  getReactionRates(T,
                   &C[0],
                   &netOut[0],
                   &createRateWorkspace[0],
                   &destroyRateWorkspace[0],
                   &stepROPWorkspace[0]);
}

void mechanism::getReactionRates_perturbROP(const double T,
                                            const double C[],
                                            const double perturbMult[],
                                            double netOut[],
                                            double createOut[],
                                            double destroyOut[],
                                            double stepOut[])
{
  perfNet->calcRatesFromTC_perturbROP(T,&C[0],&perturbMult[0],
				      &netOut[0],&createOut[0],&destroyOut[0],
				      &stepOut[0]);
}

void mechanism::getReactionRatesLimiter(const double T,
                                        const double C[],
		                        const double step_limiter[],
			                double netOut[],
                                        double createOut[],
			                double destroyOut[],
                                        double stepOut[])
{
  perfNet->calcRatesFromTC_StepLimiter(T,
                                       &C[0],
                                       &step_limiter[0],
				       &netOut[0],
                                       &createOut[0],
                                       &destroyOut[0],
				       &stepOut[0]);
}

void mechanism::getReactionRatesLimiter_perturbROP(const double T,
                                                   const double C[],
                                                   const double step_limiter[],
                                                   const double perturbMult[],
                                                   double netOut[],
                                                   double createOut[],
                                                   double destroyOut[],
                                                   double stepOut[])
{
  perfNet->calcRatesFromTC_StepLimiter_perturbROP(T,
                                                  &C[0],
                                                  &step_limiter[0],
                                                  &perturbMult[0],
                                                  &netOut[0],
                                                  &createOut[0],
                                                  &destroyOut[0],
                                                  &stepOut[0]);
}

void mechanism::getEnthalpy_RT_mr(const int nReactors, const double T[], double h_RT[]) const
{
  thermo->getH_RT_mr(nReactors,T,h_RT);
}

void mechanism::getIntEnergy_RT_mr(const int nReactors, const double T[], double u_RT[]) const
{
  thermo->getH_RT_mr(nReactors,T,u_RT);
  for(int j=0; j<nSpc*nReactors; j++)
    {u_RT[j]-=1.0;}
}


void mechanism::getMassCpFromTY_mr(const int nReactors, const double T[], const double y[],
                                      double cpSpc[], double cpReactors[]) const
{
  memset(cpReactors,0,sizeof(double)*nReactors);
  thermo->getCp_R_mr(nReactors,T,&cpSpc[0]);
  for(int j=0; j<nSpc; j++)
  {
      for(int k=0; k<nReactors; ++k)
      {
          cpSpc[nReactors*j+k]*=RuInvMolWt[j];
          cpReactors[k]+=cpSpc[nReactors*j+k]*y[nReactors*j+k];
      }
  }
}


void mechanism::getMassCvFromTY_mr(const int nReactors, const double T[], const double y[],
                                      double cvSpc[], double cvReactors[]) const
{
  memset(cvReactors,0,sizeof(double)*nReactors);
  thermo->getCp_R_mr(nReactors,T,&cvSpc[0]);
  for(int j=0; j<nSpc; j++)
  {
      for(int k=0; k<nReactors; ++k)
      {
          cvSpc[nReactors*j+k]-=1.0;  // Cv = Cp
          cvSpc[nReactors*j+k]*=RuInvMolWt[j];
          cvReactors[k]+=cvSpc[nReactors*j+k]*y[nReactors*j+k];
      }
  }
}


void mechanism::getCfromVY_mr(const int nReactors, const double v[], const double y[],
                               double c[]) const
{
  std::vector<double> dens(nReactors);
  for(int k=0; k<nReactors; ++k)
    { dens[k] = 1.0/v[k]; }
  for(int j=0; j<nSpc; j++)
  {
      for(int k=0; k<nReactors; ++k)
      {
          c[nReactors*j+k]=dens[k]*invMolWt[j]*y[nReactors*j+k];
      }
  }
}


void mechanism::getMolWtMixFromY_mr(const int nReactors, const double y[], double *mwMix) const
{
  memset(mwMix,0,sizeof(double)*nReactors);

  for(int j=0; j<nSpc; j++)
  {
    for(int k=0; k<nReactors; ++k)
      {mwMix[k]+=y[nReactors*j+k]*invMolWt[j];}
  }
  for(int k=0; k<nReactors; ++k)
      {mwMix[k]= 1.0/mwMix[k];}
}

void mechanism::getDensityFromTPY_mr(const int nReactors, const double *T, const double *P,
                                        const double y[], double *dens) const
{
  getMolWtMixFromY_mr(nReactors, y, dens);
  for(int k=0; k<nReactors; ++k)
  {
      dens[k] *= P[k]/(Ru*T[k]);
  }
}

double mechanism::getMolarAtomicOxygenRemainder(const double x[]) const
{
  int j;
  double oxygenRem=0.0;
  for(j=0; j<nSpc; j++)
    {
      oxygenRem +=x[j]*(speciesList[j].getCountOfAtomZ(ATOMIC_NUM_O)
			-2.0*speciesList[j].getCountOfAtomZ(ATOMIC_NUM_C)
			-0.5*speciesList[j].getCountOfAtomZ(ATOMIC_NUM_H));
    }
  return oxygenRem;
}

void mechanism::getMolarIdealLeanExhaust(const double xInit[],
					 double xFinal[]) const
{
  int j;
  double oxygenRem=getMolarAtomicOxygenRemainder(&xInit[0]);

  int numH,numC,numO,numN,totAtoms;
  double totH,totC,totO,totN,finalInvSum;
  int idx_co2,idx_o2,idx_n2,idx_h2o;


  if(oxygenRem < 0.0)
    {
      printf("WARNING: we do not currently have a means to estimate the\n");
      printf("         products of fuel-rich combustion. Setting final\n");
      printf("         state to initial state\n");

      for(j=0; j<nSpc; j++)
	{xFinal[j]=xInit[j];}
    }

  // find the key exhaust species indexes, and total moles of each element
  // per unit mole of the initial composition
  totH=totC=totO=totN=0.0;
  idx_o2=idx_n2=idx_h2o=idx_co2=-1;
  for(j=0; j<nSpc; j++)
    {
      totAtoms=speciesList[j].getTotalAtoms();
      numH=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_H);
      numC=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_C);
      numN=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_N);
      numO=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_O);

      totH+=xInit[j]*(double)numH;
      totC+=xInit[j]*(double)numC;
      totN+=xInit[j]*(double)numN;
      totO+=xInit[j]*(double)numO;

      if(totAtoms == 2 && numO == 2)
	{idx_o2 = j;}

      if(totAtoms == 2 && numN == 2)
	{idx_n2 = j;}

      if(totAtoms == 3 && numC == 1 && numO == 2)
	{idx_co2 = j;}

      if(totAtoms == 3 && numH == 2 && numO == 1)
	{idx_h2o = j;}
    }

  if(idx_o2 == -1 || idx_n2 == -1 || idx_co2 == -1 || idx_h2o == -1)
    {
      printf("WARNING: we are missing one of the ideal exhaust species\n");
      printf("           o2  index = %d\n",idx_o2);
      printf("           n2  index = %d\n",idx_n2);
      printf("           co2 index = %d\n",idx_co2);
      printf("           h2o index = %d\n",idx_h2o);
      printf("         Setting final composition to initial composition\n");

      for(j=0; j<nSpc; j++)
	{xFinal[j]=xInit[j];}
    }
  else
    {
      for(j=0; j<nSpc; j++)
	{xFinal[j]=0.0;}

      // the following assignments of complete combustion are based on the
      // assumption of there existing excess oxygen
      xFinal[idx_n2]  = 0.5*totN;      // all atomic nitrogen goes to N2
      xFinal[idx_h2o] = 0.5*totH;      // all atomic hydrogen goes to H2
      xFinal[idx_co2] = totC;          // all atomic carbon goes to CO2
      xFinal[idx_o2]  = 0.5*oxygenRem; // all remaining oxygen goes to O2

      // current composition is the number of moles of a species per mole of
      // initial composition, so renormalize to ideal exhaust composition
      finalInvSum=1.0/(xFinal[idx_n2]+xFinal[idx_h2o]+
		       xFinal[idx_co2]+xFinal[idx_o2]);

      xFinal[idx_n2]  *= finalInvSum;
      xFinal[idx_h2o] *= finalInvSum;
      xFinal[idx_co2] *= finalInvSum;
      xFinal[idx_o2]  *= finalInvSum;
    }
}
void mechanism::getMolarIdealExhaust(const double xInit[],
					 double xFinal[]) const
{
  int j;
  double oxygenRem=getMolarAtomicOxygenRemainder(&xInit[0]);

  int numH,numC,numO,numN,totAtoms;
  double totH,totC,totO,totN,finalInvSum;
  int idx_co2,idx_o2,idx_n2,idx_h2o;
  int idx_co;

  // find the key exhaust species indexes, and total moles of each element
  // per unit mole of the initial composition
  totH=totC=totO=totN=0.0;
  idx_o2=idx_n2=idx_h2o=idx_co2=idx_co=-1;
  for(j=0; j<nSpc; j++)
    {
      totAtoms=speciesList[j].getTotalAtoms();
      numH=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_H);
      numC=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_C);
      numN=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_N);
      numO=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_O);

      totH+=xInit[j]*(double)numH;
      totC+=xInit[j]*(double)numC;
      totN+=xInit[j]*(double)numN;
      totO+=xInit[j]*(double)numO;

      if(totAtoms == 2 && numO == 2)
	{idx_o2 = j;}

      if(totAtoms == 2 && numN == 2)
	{idx_n2 = j;}

      if(totAtoms == 3 && numC == 1 && numO == 2)
	{idx_co2 = j;}

      if(totAtoms == 2 && numC == 1 && numO == 1)
	{idx_co  = j;}

      if(totAtoms == 3 && numH == 2 && numO == 1)
	{idx_h2o = j;}
    }

  if(idx_o2 == -1 || idx_n2 == -1 || idx_co2 == -1 || idx_h2o == -1 ||
     idx_co == -1)
    {
      printf("WARNING: we are missing one of the ideal exhaust species\n");
      printf("           o2  index = %d\n",idx_o2);
      printf("           n2  index = %d\n",idx_n2);
      printf("           co  index = %d\n",idx_co);
      printf("           co2 index = %d\n",idx_co2);
      printf("           h2o index = %d\n",idx_h2o);
      printf("         Setting final composition to initial composition\n");

      for(j=0; j<nSpc; j++)
	{xFinal[j]=xInit[j];}
    }
  else
    {
      for(j=0; j<nSpc; j++)
	{xFinal[j]=0.0;}

      if(oxygenRem >= 0.0)
	{
	  // the following assignments of complete combustion are based on the
	  // assumption of there existing excess oxygen
	  xFinal[idx_n2]  = 0.5*totN;      // all atomic nitrogen goes to N2
	  xFinal[idx_h2o] = 0.5*totH;      // all atomic hydrogen goes to H2
	  xFinal[idx_co2] = totC;          // all atomic carbon goes to CO2
	  xFinal[idx_o2]  = 0.5*oxygenRem; // all remaining oxygen goes to O2

	  // current composition is the number of moles of a species per mole of
	  // initial composition, so renormalize to ideal exhaust composition
	  finalInvSum=1.0/(xFinal[idx_n2]+xFinal[idx_h2o]+
			   xFinal[idx_co2]+xFinal[idx_o2]);

	  xFinal[idx_n2]  *= finalInvSum;
	  xFinal[idx_h2o] *= finalInvSum;
	  xFinal[idx_co2] *= finalInvSum;
	  xFinal[idx_o2]  *= finalInvSum;
	}
      else
	{
	  // recalculate oxygen remainder assuming H -> H2O and C -> CO
	  oxygenRem = totO - totC - 0.5*totH;
	  if(oxygenRem >= 0.0)
	    {
	      xFinal[idx_n2]  = 0.5*totN;    // all atomic nitrogen goes to N2
	      xFinal[idx_h2o] = 0.5*totH;    // all atomic hydrogen goes to H2
	      xFinal[idx_co]  = totC-oxygenRem;
	      xFinal[idx_co2] = oxygenRem;

	      // current composition is the number of moles of a species per
	      // mole of initial composition, so renormalize to ideal exhaust
	      // composition
	      finalInvSum=1.0/(xFinal[idx_n2]+xFinal[idx_h2o]+
			       xFinal[idx_co2]+xFinal[idx_co]);

	      xFinal[idx_n2]  *= finalInvSum;
	      xFinal[idx_h2o] *= finalInvSum;
	      xFinal[idx_co2] *= finalInvSum;
	      xFinal[idx_co]  *= finalInvSum;
	    }
	  else
	    {
	      printf("WARNING: we do not currently have an ideal exhaust\n");
	      printf("         model when there is not enough oxygen to\n");
	      printf("         reduce all the C to CO.\n");
	      printf("         Setting final composition to initial composition\n");
	      for(j=0; j<nSpc; j++)
		{xFinal[j]=xInit[j];}
	    }
	}
    }
}

double mechanism::getProgressEquivalenceRatio(const double y[]) const
{
  const double z_prime = 0.0; //percent of oxygen in fuel
                              //if important needs to be set per mech
  double totH=0.0;
  double totC=0.0;
  double totO=0.0;
  for(int j=0; j<nSpc; j++)
  {
      int totAtoms=speciesList[j].getTotalAtoms();
      int numH=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_H);
      int numC=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_C);
      int numO=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_O);

      if( !(totAtoms == 3 && numC == 1 && numO == 2) && //CO2
          !(totAtoms == 3 && numH == 2 && numO == 1) )  //H2O
      {
        totH+=y[j]*(double)numH*invMolWt[j];
        totC+=y[j]*(double)numC*invMolWt[j];
        totO+=y[j]*(double)numO*invMolWt[j];
      }
  }
  return ((2.0-z_prime)*totC  + 0.5*totH)/(totO-z_prime*totC);
}


void mechanism::getModifiedEquivalenceRatios(const double y[],
                                             double &smallmxf,
                                             double &bigmxf) const
{
  double totH=0.0;
  double totC=0.0;
  double totO=0.0;
  double totH2=0.0;
  double totC2=0.0;
  double totO2=0.0;
  for(int j=0; j<nSpc; j++)
  {
      int totAtoms=speciesList[j].getTotalAtoms();
      int numH=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_H);
      int numC=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_C);
      int numO=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_O);

      totH+=y[j]*(double)numH*invMolWt[j];
      totC+=y[j]*(double)numC*invMolWt[j];
      totO+=y[j]*(double)numO*invMolWt[j];
      if( !(totAtoms == 3 && numC == 1 && numO == 2) && //CO2
          !(totAtoms == 3 && numH == 2 && numO == 1) )  //H2O
      {
        totH2+=y[j]*(double)numH*invMolWt[j];
        totC2+=y[j]*(double)numC*invMolWt[j];
        totO2+=y[j]*(double)numO*invMolWt[j];
      }
  }
  double chfact = 2.0*totC + 0.5*totH;
  double xu;
  if( totO >= chfact ) //lean
  {
    bigmxf = 0.5*chfact/totO;
    xu = chfact > 0.0 ? (2.0*totC2+0.5*totH2)/chfact : 0.0;
  } else {            //rich
    bigmxf = 1.0 - 0.5*totO/chfact;
    xu = totO > 0.0 ? totO2/totO : 0.0;
  }
  smallmxf = bigmxf * xu;
}

double mechanism::getCHValue(const double y[]) const
{
  double totH=0.0;
  double totC=0.0;
  for(int j=0; j<nSpc; j++)
  {
      int totAtoms=speciesList[j].getTotalAtoms();
      int numH=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_H);
      int numC=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_C);
      int numO=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_O);

      if( !(totAtoms == 3 && numC == 1 && numO == 2) && //CO2
          !(totAtoms == 3 && numH == 2 && numO == 1) )  //H2O
      {
        totH+=y[j]*(double)numH*invMolWt[j];
        totC+=y[j]*(double)numC*invMolWt[j];
      }
  }
  return (2.0*totC  + 0.5*totH);
}

void mechanism::getOxygenAtomCount(int numO[]) const
{
  for(int j=0; j<nSpc; ++j)
  {
    numO[j]=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_O);
  }
}

void mechanism::getCarbonAtomCount(int numC[]) const
{
  for(int j=0; j<nSpc; ++j)
  {
    numC[j]=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_C);
  }
}

void mechanism::getHydrogenAtomCount(int numH[]) const
{
  for(int j=0; j<nSpc; ++j)
  {
    numH[j]=speciesList[j].getCountOfAtomZ(ATOMIC_NUM_H);
  }
}

void mechanism::getSpeciesHydrogenCount(int num_atoms[]) const
{
  const int num_species = nSpc;
  for(int j=0; j<num_species; ++j) {
    num_atoms[j] = speciesList[j].getCountOfAtomZ(ATOMIC_NUM_H);
  }
}

void mechanism::getSpeciesNitrogenCount(int num_atoms[]) const
{
  const int num_species = nSpc;
  for(int j=0; j<num_species; ++j) {
    num_atoms[j] = speciesList[j].getCountOfAtomZ(ATOMIC_NUM_N);
  }
}
void mechanism::getSpeciesCarbonCount(int num_atoms[]) const
{
  const int num_species = nSpc;
  for(int j=0; j<num_species; ++j) {
    num_atoms[j] = speciesList[j].getCountOfAtomZ(ATOMIC_NUM_C);
  }
}
void mechanism::getSpeciesOxygenCount(int num_atoms[]) const
{
  const int num_species = nSpc;
  for(int j=0; j<num_species; ++j) {
    num_atoms[j] = speciesList[j].getCountOfAtomZ(ATOMIC_NUM_O);
  }
}

void mechanism::getSpeciesArgonCount(int num_atoms[]) const
{
  const int num_species = nSpc;
  for(int j=0; j<num_species; ++j) {
    num_atoms[j] = speciesList[j].getCountOfAtomZ(ATOMIC_NUM_AR);
  }
}

void mechanism::getSpeciesHeliumCount(int num_atoms[]) const
{
  const int num_species = nSpc;
  for(int j=0; j<num_species; ++j) {
    num_atoms[j] = speciesList[j].getCountOfAtomZ(ATOMIC_NUM_HE);
  }
}

void mechanism::initExternalFuncs()
{
#ifndef _WIN32
  char * externalFuncLibPath;
  char * dl_error;
  char * ext_flags;
  unsigned char ext_active_flags = 0; //Make user turn on explicitly
  externalFuncLibHandle = NULL;

  ext_flags = getenv("ZERORK_EXT_FLAGS");
  if (ext_flags!=NULL)
  {
     ext_active_flags = atoi(ext_flags);
  }
  if( ext_active_flags <= 0 ) return;

  externalFuncLibPath = getenv("ZERORK_EXT_FUNC_LIB");
  if (externalFuncLibPath!=NULL)
  {
      printf("ZERORK_MECHANISM: $ZERORK_EXT_FUNC_LIB is: %s\n",externalFuncLibPath);
      externalFuncLibHandle = dlopen(externalFuncLibPath,RTLD_NOW);
      if (!externalFuncLibHandle)
      {
         fprintf(stderr, "%s\n", dlerror());
         return;
      }
  } else {
    return;
  }

  if( externalFuncLibHandle != NULL ) // Attempt to load from external library.
  {
      ex_func_check = (external_func_check_t) dlsym(externalFuncLibHandle, "external_func_check");
      if ((dl_error = dlerror()) != NULL)
      {
         fprintf(stderr, "%s\n", dl_error);
         return;
      }
      if(!ex_func_check(nSpc,nStep)) {
        return;
      }

      ex_func_calc_rates = (external_func_rates_t) dlsym(externalFuncLibHandle, "external_func_rates");
      if ((dl_error = dlerror()) != NULL)
      {
         fprintf(stderr, "%s\n", dl_error);
         return;
      }

      ex_func_calc_arrh = (external_func_arrh_t) dlsym(externalFuncLibHandle, "external_func_arrh");
      if ((dl_error = dlerror()) != NULL)
      {
         fprintf(stderr, "%s\n", dl_error);
         return;
      }

      ex_func_calc_keq = (external_func_keq_t) dlsym(externalFuncLibHandle, "external_func_keq");
      if ((dl_error = dlerror()) != NULL)
      {
         fprintf(stderr, "%s\n", dl_error);
         return;
      }
  } else {
    return;
  }

  if( ext_active_flags & 0x01 )
  {
      printf("ZERORK_MECHANISM: activating external rates.\n");
      perfNet->setUseExRates();
      perfNet->setExRatesFunc(ex_func_calc_rates);
  }
  if( ext_active_flags & 0x02 )
  {
      printf("ZERORK_MECHANISM: activating external arrh.\n");
      Kconst->setUseExArrh();
      Kconst->setExArrhFunc(ex_func_calc_arrh);
  }
  if( ext_active_flags & 0x04 )
  {
      printf("ZERORK_MECHANISM: activating external keq.\n");
      Kconst->setUseExKeq();
      Kconst->setExKeqFunc(ex_func_calc_keq);
  }
#endif
}

void mechanism::buildReactionString(const int idx,
                                    string &str)
{
  int stepId,falloffFlag;
  string spcStr;

  assert(infoNet != NULL);
  assert(speciesList != NULL);

  stepId = getStepIdxOfRxn(idx,1);
  str.clear();

  // store the first reactant species
  spcStr = getSpeciesName(getSpecIdxOfStepReactant(stepId,0));
  str = spcStr;
  // store the remaining reactant species
  for(int j=1; j<infoNet->getOrderOfStep(stepId); j++) {
    spcStr = getSpeciesName(getSpecIdxOfStepReactant(stepId,j));
    str += " + " + spcStr;
  }
  // check for third body
  if(infoNet->getThirdBodyFlagOfReaction(idx) == 1) {
    str += " + M";
  }
  // check for falloff
  falloffFlag = infoNet->getFalloffFlagOfReaction(idx);
  if(falloffFlag == -1) {
    str += " (+M)";
  }
  if(falloffFlag > 0) {
    spcStr = getSpeciesName(falloffFlag-1);
    str += " (+" + spcStr + ")";
  }
  // check for reversible reactions
  if(infoNet->getReversibleFlagOfReaction(idx) == 1) {
    str += " <=> ";
  }
  else {
    str += " => ";
  }

  // store the first product species
  spcStr = getSpeciesName(getSpecIdxOfStepProduct(stepId,0));
  str += spcStr;

  // store the remaining product species
  for(int j=1; j<infoNet->getNumProductsOfStep(stepId); j++) {
    spcStr = getSpeciesName(getSpecIdxOfStepProduct(stepId,j));
    str += " + " + spcStr;
  }
  // check for third body
  if(infoNet->getThirdBodyFlagOfReaction(idx) == 1) {
    str += " + M";
  }
  // check for falloff
  falloffFlag = infoNet->getFalloffFlagOfReaction(idx);
  if(falloffFlag == -1) {
    str += " (+M)";
  }
  if(falloffFlag > 0) {
    spcStr = getSpeciesName(falloffFlag-1);
    str += " (+" + spcStr + ")";
  }
}

void mechanism::generateCircosFilesTPY(const double T, const double P,
                                       const double y[],
                                       const int num_reacs,
                                       const int exclude_unimolecular,
                                       const char *filename_prefix
                                       )
{
   double specific_volume = 1.0/getDensityFromTPY(T, P, y);
   std::vector<double> concentration(nSpc);
   std::vector<double> netROP(nSpc);
   getCfromVY(specific_volume, y, &concentration[0]);
   getReactionRates(T, &concentration[0],
                    &netROP[0],
                    &createRateWorkspace[0],
                    &destroyRateWorkspace[0],
                    &stepROPWorkspace[0]);


    int numReactants, numProducts;
    double maxROP,minROP;
    maxROP = 0.0;
    minROP = 1.0e+300;
    //Find unscaled max and min
    for(int j=0; j<nStep; ++j)
    {
        numReactants = getOrderOfStep(j);
        numProducts = getNumProductsOfStep(j);
        if(exclude_unimolecular &&
           (numReactants == 1 || numProducts == 1)) continue;
        if(stepROPWorkspace[j] > maxROP) maxROP=stepROPWorkspace[j];
        if(stepROPWorkspace[j] < minROP) minROP=stepROPWorkspace[j];
    }

    //Scale by max and take log10
    std::vector<double> scaledROPS(nStep);
    for(int j=0; j<nStep; ++j)
    {
        numReactants = getOrderOfStep(j);
        numProducts = getNumProductsOfStep(j);
        if(exclude_unimolecular &&
           (numReactants == 1 || numProducts == 1))
        {
            scaledROPS[j] = -1.0e+30;
            continue;
        }
        scaledROPS[j] = stepROPWorkspace[j]/maxROP;
        if(scaledROPS[j] >= 0.0)  //negative concentrations can yeild negative ROPs
        {
            scaledROPS[j] = log10(scaledROPS[j]);
        }
        else
        {
            scaledROPS[j] = -1.0e+30;
        }
    }

    maxROP = -1.0e+300;
    minROP = 1.0e+300;
    //Find scaled max and min
    for(int j=0; j<nStep; ++j)
    {
        numReactants = getOrderOfStep(j);
        numProducts = getNumProductsOfStep(j);
        if(numReactants == 1 || numProducts == 1) continue;
        if(scaledROPS[j] > maxROP) maxROP=scaledROPS[j];
        if(scaledROPS[j] < minROP) minROP=scaledROPS[j];
    }

    std::vector<double> sortedROPS(nStep);

    memcpy(&sortedROPS[0],&scaledROPS[0],sizeof(double)*nStep);
    std::sort(sortedROPS.begin(),sortedROPS.end());

    int nkeep=min(num_reacs,nStep);
    minROP=sortedROPS[nkeep-1];

    int iBand = 0;
    int total_bands = getTotalProducts()+getTotalReactants();
    std::vector<int> band_id(total_bands);
    std::vector<int> band_spc(total_bands);
    std::vector<int> band_start(total_bands);
    std::vector<int> band_end(total_bands);
    std::vector<int> band_dir(total_bands);
    for(int i=0; i< nSpc; ++i)
    {
        int nSpecieRxns = 0;
        for (int j=0; j<getNumSteps(); ++j)
        {
            if(minROP > scaledROPS[j]) continue;
            numReactants = getOrderOfStep(j);
            numProducts = getNumProductsOfStep(j);
            if(exclude_unimolecular &&
               (numReactants == 1 || numProducts == 1)) continue;
            for(int k=0; k<numReactants; ++k)
            {
                if(getSpecIdxOfStepReactant(j,k) == i)
                {
                    band_id[iBand] = j;
                    band_spc[iBand] = i;
                    band_start[iBand] = nSpecieRxns;
                    band_end[iBand] = nSpecieRxns+1;
                    band_dir[iBand] = +1;
                    ++iBand;
                    ++nSpecieRxns;
                    break; //only count species once
                }
            }
            for(int k=0; k<numProducts; ++k)
            {
                if(getSpecIdxOfStepProduct(j,k) == i)
                {
                    band_id[iBand] = j;
                    band_spc[iBand] = i;
                    band_start[iBand] = nSpecieRxns;
                    band_end[iBand] = nSpecieRxns+1;
                    band_dir[iBand] = -1;
                    ++nSpecieRxns;
                    ++iBand;
                    break; //only count species once
                }
            }
        }
    }


    int iLink = 0;
    int total_links = 2*total_bands;
    std::vector<int> link_rxn(total_links);
    std::vector<int> link_reac(total_links);
    std::vector<int> link_prod(total_links);
    std::vector<int> link_dir(total_links);
    std::vector<int> link_rxnNum(total_links);
    std::vector<int> specieRxnCounts(nSpc,0);

    for (int j=0; j < getNumSteps(); ++j)
    {
        if(minROP > scaledROPS[j]) continue; //excluded by rate
        numReactants = getOrderOfStep(j);
        numProducts = getNumProductsOfStep(j);
        if(exclude_unimolecular
           && (numReactants == 1 || numProducts == 1)) continue;
        int lastReactant = -1;
        for(int iR = 0; iR<numReactants; ++iR)
        {
            int currReac = getSpecIdxOfStepReactant(j,iR);
            if( currReac == lastReactant ) continue;
            specieRxnCounts[currReac] += 1;
            lastReactant = currReac;
        }
        int lastProduct = -1;
        for(int iP = 0; iP<numProducts; ++iP)
        {
            int currProd = getSpecIdxOfStepProduct(j,iP);
            if( currProd == lastProduct ) continue;
            specieRxnCounts[currProd] += 1;
            lastProduct = currProd;
        }
        lastReactant = -1;
        for(int iR = 0; iR<numReactants; ++iR)
        {
            int currReac = getSpecIdxOfStepReactant(j,iR);
            if( currReac == lastReactant ) continue; //don't double count
            lastReactant = currReac;
            lastProduct = -1;
            for(int iP = 0; iP<numProducts; ++iP)
            {
                int currProd = getSpecIdxOfStepProduct(j,iP);
                if( currProd == lastProduct ) continue; //don't double count;
                lastProduct = currProd;

                link_rxn[iLink] = j;
                link_reac[iLink] = currReac;
                link_prod[iLink] = currProd;
                link_dir[iLink] = +1;
                link_rxnNum[iLink] = specieRxnCounts[currReac]-1;

                iLink++;

                link_rxn[iLink] = j;
                link_reac[iLink] = currReac;
                link_prod[iLink] = currProd;
                link_dir[iLink] = -1;
                link_rxnNum[iLink] = specieRxnCounts[currProd]-1;

                iLink++;
            }
        }
    }


    char karyotype_filename[64];
    snprintf(karyotype_filename,64,"%s.karyo.txt",filename_prefix);
    FILE * karyotype_file = fopen(karyotype_filename,"w");
    //Karyotypes
    for(int i=0; i< nSpc; ++i)
    {
        if(specieRxnCounts[i] == 0) continue;
        fprintf(karyotype_file,"chr - %s %s 0 %d green\n",
                getSpeciesName(i),getSpeciesName(i),specieRxnCounts[i]);
    }

    //Bands
    for(int j = 0; j < iBand; ++j)
    {
        if(band_dir[j] > 0) //destroy
        {
            fprintf(karyotype_file,"band %s reac%d reac%d %d %d red\n",
                    getSpeciesName(band_spc[j]),band_id[j],band_id[j],
                    band_start[j],band_end[j]);
        }
        else //create
        {
            fprintf(karyotype_file,"band %s prod%d prod%d %d %d green\n",
                    getSpeciesName(band_spc[j]),band_id[j],band_id[j],
                    band_start[j],band_end[j]);
        }
    }
    fclose(karyotype_file);

    char links_filename[64];
    snprintf(links_filename,64,"%s.links.txt",filename_prefix);
    FILE * links_file = fopen(links_filename,"w");
    //Links
    for(int j = 0; j < iLink; ++j)
    {
        int hueIdx = 0;
        if(scaledROPS[link_rxn[j]] < minROP){ hueIdx = 240; } //blue
        else if(scaledROPS[link_rxn[j]] > maxROP){ hueIdx = 0;} //red
        else { hueIdx = 240 - int(240*(scaledROPS[link_rxn[j]]-minROP)/(maxROP-minROP)); }
        if(link_dir[j] > 0)
        {
           fprintf(links_file,"rxn%dreac%dprod%d %s %d %d color=hue%03d_a2,z=-%d\n",
                   link_rxn[j],link_reac[j],link_prod[j],
                   getSpeciesName(link_reac[j]),
                   link_rxnNum[j],link_rxnNum[j]+1,hueIdx,hueIdx);
        }
        else
        {
           fprintf(links_file,"rxn%dreac%dprod%d %s %d %d color=hue%03d_a2,z=-%d\n",
                   link_rxn[j],link_reac[j],link_prod[j],
                   getSpeciesName(link_prod[j]),
                   link_rxnNum[j],link_rxnNum[j]+1,hueIdx,hueIdx);
        }
    }
    fclose(links_file);
}

} // namespace zerork
