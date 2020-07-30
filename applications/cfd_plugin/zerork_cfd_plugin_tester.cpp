
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utility_funcs.h"
#include "ZeroRKCFDPluginTesterIFP.h"

#include "zerork/mechanism.h"
#include "zerork_cfd_plugin.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

template <typename T>
static std::vector<T> linspace(T a, T b, size_t N);

static void getFracsFromCompMap(zerork::mechanism& mech,
    const std::map<std::string, double> comp, double fracArray[]);

static void setReactorMassFrac(int nReactors, const zerork::mechanism& mech,
    double systemState[], const double fuelMoleFrac[],
    const double oxidMoleFrac[], const double reactorPhi[],
    const double reactorEgr[]);


void zerork_cfd_plugin_tester(int inp_argc, char **inp_argv)
{
  if(inp_argc != 2)
    {
      printf("ERROR: incorrect command line usage\n");
      printf("       use instead %s <idt sweep input>\n",inp_argv[0]);
      fflush(stdout); exit(-1);
    }

  ZeroRKCFDPluginTesterIFP inputFileDB(inp_argv[1]);

  int nReactors = inputFileDB.nReactors();

  zerork::mechanism mech(
      inputFileDB.mechFile().c_str(),
      inputFileDB.thermFile().c_str(),
      inputFileDB.mechLogFile().c_str(),
      1 // verbosity
  );

  int nSpc=mech.getNumSpecies();
  int nState=nSpc+1;
  int nStep=mech.getNumSteps();

  double cv_maxsteps = inputFileDB.maxSteps();
  double cv_maxord   = inputFileDB.maxOrd();
  double cv_rtol     = inputFileDB.relTol();
  double cv_atol     = inputFileDB.absTol();
  int abstol_dens = 0;
  double cv_maxdt_internal = inputFileDB.maxDtInternal();
  double cv_sparsethresh = inputFileDB.precThresh();
  double cv_nlconvcoef = inputFileDB.cvNlConvCoeff();
  double cv_epslin = inputFileDB.cvEpsLin();
  const char* cklogfile = inputFileDB.mechLogFile().c_str();
  const char* reactorlogfile = inputFileDB.outFile().c_str();
  int multireac = 0;
  int gpu_id = -1;

  int lsolver = -1;
  if( inputFileDB.linearSolver() == std::string("DirectDense") ) {
    lsolver = 0;
  } else if( inputFileDB.linearSolver() == std::string("DirectDenseDVD") ) {
    lsolver = 1;
  } else if( inputFileDB.linearSolver() == std::string("IterativeSparse") ) {
    lsolver = 2;
  }

  int verbosity = 5;

  //Setup reactor_lib
  zerork_cfd_plugin_setup_full(verbosity,cv_maxsteps,cv_maxord,cv_rtol,cv_atol,
                         abstol_dens,
                         cv_maxdt_internal,cv_sparsethresh,cv_nlconvcoef,
                         cv_epslin,lsolver,
                         inputFileDB.mechFile().c_str(),
                         inputFileDB.thermFile().c_str(),
                         cklogfile,reactorlogfile,
                         &multireac, &gpu_id);


  //Set up reactor initial states
  std::vector<double> reactorT(nReactors);
  std::vector<double> reactorP(nReactors);
  std::vector<double> reactorMassFrac(nReactors*nSpc);
  std::vector<double> reactorCost(nReactors);
  std::vector<double> reactorGpu(nReactors);

  if(inputFileDB.stateFiles().size() != 0)
  {
    //Load initial state from files
    if(inputFileDB.stateFiles().size() != nReactors)
    {
      printf("stateFiles present but not matching nReactors. Quitting.\n");
      exit(-1);
    }
    for(int k=0; k<nReactors; ++k)
    {
      FILE* stateFile = fopen(inputFileDB.stateFiles()[k].c_str(),"r");
      if(stateFile == NULL)
      {
        printf("Unable to read file: %s\n",inputFileDB.stateFiles()[k].c_str());
        exit(-1);
      }
      double val;
      for(int j=0;j<nSpc;++j)
      {
        int matches = fscanf(stateFile,"%lf",&val);
        if( matches != 1 )
        {
          printf("Failed to find value in file %s\n",
                 inputFileDB.stateFiles()[k].c_str());
          exit(-1);
        }
        reactorMassFrac[k*nSpc+j] = val;
      }
      fscanf(stateFile,"%lf",&val);
      reactorT[k] = val;
      fscanf(stateFile,"%lf",&val);
      reactorP[k] = val;
      fclose(stateFile);
    }
  }
  else
  {
     //Set initial states using linear ranges of T,P,phi, and egr
    double phiMin, phiMax, egrMin, egrMax, TMin, TMax, pMin, pMax;
    phiMin = inputFileDB.phiMin();
    phiMax = inputFileDB.phiMax();
    egrMin = inputFileDB.egrMin();
    egrMax = inputFileDB.egrMax();
    TMin = inputFileDB.TMin();
    TMax = inputFileDB.TMax();
    pMin = inputFileDB.pMin();
    pMax = inputFileDB.pMax();

    reactorT = linspace(TMin,TMax,nReactors);
    reactorP = linspace(pMin,pMax,nReactors);
    std::vector<double> reactorPhi = linspace(phiMin,phiMax,nReactors);
    std::vector<double> reactorEgr = linspace(egrMin,egrMax,nReactors);

    std::vector<double> initFuelMassFrac(nSpc);
    std::vector<double> initOxidMassFrac(nSpc);
    getFracsFromCompMap(mech,inputFileDB.fuelComp(),&initFuelMassFrac[0]);
    getFracsFromCompMap(mech,inputFileDB.oxidizerComp(),&initOxidMassFrac[0]);

    setReactorMassFrac(nReactors, mech, &reactorMassFrac[0],
                       &initFuelMassFrac[0], &initOxidMassFrac[0],
                       &reactorPhi[0], &reactorEgr[0]);
  }


  double tend = inputFileDB.reactorTime();

  // timing data
  double startTime,stopTime,otherTime,simTime;
  startTime=getHighResolutionTime();
  zerork_reset_stats();
  zerork_solve_reactors(nReactors, tend, &reactorT[0], &reactorP[0],
                     &reactorMassFrac[0], &reactorCost[0], &reactorGpu[0]);
  zerork_print_stats();

  stopTime=getHighResolutionTime();
  simTime=stopTime-startTime;
  printf("simTime : %g s\n",simTime);
}

int main(int argc, char **argv)
{
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
  zerork_cfd_plugin_tester(argc, argv);
#ifdef USE_MPI
  MPI_Finalize();
#endif
  exit(0);
}


template <typename T>
static std::vector<T> linspace(T a, T b, size_t N)
{
  std::vector<T> xs(N);
  T h = (b - a) / static_cast<T>(N-1);
  T val = a;
  typename std::vector<T>::iterator x;
  for(x = xs.begin(); x != xs.end(); ++x)
  {
    *x = val;
    val += h;
  }
  return xs;
}

template <typename T>
static T normalize(const size_t N, T v[])
{
  size_t j;
  T sum=0.0;
  T invSum=1.0;

  for(j=0; j<N; ++j)
    {sum+=v[j];}

  if(sum != 0.0)
    {invSum=1.0/sum;}

  for(j=0; j<N; ++j)
    {v[j]*=invSum;}
  return sum;
}

static void getFracsFromCompMap(zerork::mechanism& mech,
    const std::map<std::string, double> comp,
    double fracArray[]
)
{
  int nSpc = mech.getNumSpecies();
  memset(fracArray,0,nSpc*sizeof(double));
  std::map<std::string, double>::const_iterator mapit;
  for(mapit = comp.begin();
         mapit != comp.end(); ++mapit)
  {
    std::string spcName = mapit->first;
    int spcIdx = mech.getIdxFromName(spcName.c_str());
    if(spcIdx == -1)
    {
      std::cerr << "Species \"" << spcName << "\" in composition array \""
                << "\" not found in mechanism." << std::endl;
      exit(-1);
    }
    fracArray[spcIdx]=mapit->second;
  }
  normalize(nSpc,fracArray);
}

static void setReactorMassFrac(int nReactors, const zerork::mechanism& mech,
    double massFracs[], const double fuelMoleFrac[],
    const double oxidMoleFrac[], const double reactorPhi[],
    const double reactorEgr[])
{
  int j,k;
  int nSpc = mech.getNumSpecies();
  std::vector<double> frMole(nSpc);
  std::vector<double> frMass(nSpc);
  std::vector<double> exMass(nSpc);
  std::vector<double> exMole(nSpc);

  std::vector<double> freshMoleFrac(nReactors*nSpc);

  // at this point the fuel and oxidizer mole fraction compositions should
  // be properly normalized
  double fuelOxygenBal = mech.getMolarAtomicOxygenRemainder(&fuelMoleFrac[0]);
  double oxidOxygenBal = mech.getMolarAtomicOxygenRemainder(&oxidMoleFrac[0]);
  double moleOxidStoicPerFuel = fabs(fuelOxygenBal)/oxidOxygenBal;

  // compute the fresh mixture mole fraction based on phi and the number
  // of moles of oxidizer needed for stoichiometry per mole of fuel
  for(k=0; k<nReactors; ++k)
  {
      for(j=0; j<nSpc; j++)
      {
          freshMoleFrac[nSpc*k+j] = reactorPhi[k]*fuelMoleFrac[j]
                            +moleOxidStoicPerFuel*oxidMoleFrac[j];
          frMole[j] = freshMoleFrac[nSpc*k+j];
      }

      normalize(nSpc,&frMole[0]);
      mech.getYfromX(&frMole[0],&frMass[0]);

      if(reactorEgr[k] > 0.0)
      {
          // get the ideal exhaust molar composition
          mech.getMolarIdealExhaust(&frMole[0],&exMole[0]);

          // convert the exhaust mole fractions to mass fractions
          mech.getYfromX(&exMole[0], &exMass[0]);

          // create the initial mass fraction of the blended intake composition
          // here egr represents the fraction by mass the ideal exhaust composition
          // is in the intake composition
      }

      for(j=0; j<nSpc; j++)
      {
          massFracs[nSpc*k+j]=(1.0-reactorEgr[k])*frMass[j]+reactorEgr[k]*exMass[j];
      }
   }
}


