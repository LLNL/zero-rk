
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iomanip>
#include <sstream>

#include "utility_funcs.h"
#include "ZeroRKCFDPluginTesterIFP.h"

#include "zerork/mechanism.h"
#include "file_utilities.h"
#include "zerork_cfd_plugin.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

using zerork::getHighResolutionTime;

template <typename T>
static std::vector<T> linspace(T a, T b, size_t N);

static void getFracsFromCompMap(zerork::mechanism& mech,
    const std::map<std::string, double> comp, double fracArray[]);

static void setReactorMassFrac(int nReactors, const zerork::mechanism& mech,
    double systemState[], const double fuelMoleFrac[],
    const double oxidMoleFrac[], const double reactorPhi[],
    const double reactorEgr[]);

static void log_output(int step, double time, int n_print_reactors,
                       std::vector<double> reactorT, std::vector<double> reactorP,
                       std::vector<double> reactorDPDT, std::vector<double> reactorMassFrac,
                       std::vector<double> reactorCost, std::vector<double> reactorGpu,
                       std::vector<int> log_species_indexes, std::vector<std::string> log_species_names,
                       int nsp, std::vector<std::shared_ptr<std::ofstream>> reactor_log_files);


typedef struct UserData {
  int nsteps;
  double time;
} user_data_t;

static zerork_handle zrm_handle;
static user_data_t ud;

static int callback(int reactor_id, int nsteps, double time, double dt, const double* y, const double* ydot, void* user_data) {
  user_data_t* ud = static_cast<user_data_t*>(user_data);
  //printf("%d, %g: %g\n",nsteps, time, y[0]);
  ud->nsteps += 1;
  ud->time += dt;
  return 0;
}

void zerork_reactor(int inp_argc, char **inp_argv)
{
  if(inp_argc != 2) {
    printf("ERROR: incorrect command line usage\n");
    printf("       use instead %s <idt sweep input>\n",inp_argv[0]);
    fflush(stdout); exit(-1);
  }


  ZeroRKCFDPluginTesterIFP inputFileDB(inp_argv[1]);

  const int constant_volume = inputFileDB.constant_volume() ? 1 : 0;
  const int stop_after_ignition = inputFileDB.stop_after_ignition();
  const double delta_temperature_ignition = inputFileDB.delta_temperature_ignition();

  const char* mechfilename = inputFileDB.mechanism_file().c_str();
  const char* thermfilename = inputFileDB.thermo_file().c_str();
  int error_state = 0;
  zerork_status_t status_mech = ZERORK_STATUS_SUCCESS;
  zerork_status_t status_options = ZERORK_STATUS_SUCCESS;
  zerork_status_t status_other = ZERORK_STATUS_SUCCESS;
#ifdef USE_OMP
  #pragma omp threadprivate(zrm_handle, ud)
  #pragma omp parallel reduction(+:error_state)
  {
#endif
  zrm_handle = zerork_reactor_init();
  ud.nsteps = 0;
  ud.time = 0.0;
  if(inputFileDB.zerork_cfd_plugin_input().size() > 0) {
    status_options = zerork_reactor_read_options_file(inputFileDB.zerork_cfd_plugin_input().c_str(), zrm_handle);
  }
  status_other = zerork_reactor_set_mechanism_files(mechfilename, thermfilename, zrm_handle);
  if(status_other != ZERORK_STATUS_SUCCESS) error_state += 1;
  status_mech = zerork_reactor_load_mechanism(zrm_handle);
  status_other = zerork_reactor_set_callback_fn(callback, &ud, zrm_handle);
  if(status_other != ZERORK_STATUS_SUCCESS) error_state += 1;
  status_other = zerork_reactor_set_int_option("constant_volume", constant_volume, zrm_handle);
  if(status_other != ZERORK_STATUS_SUCCESS) error_state += 1;
  status_other = zerork_reactor_set_int_option("stop_after_ignition", stop_after_ignition, zrm_handle);
  if(status_other != ZERORK_STATUS_SUCCESS) error_state += 1;
  status_other = zerork_reactor_set_double_option("delta_temperature_ignition", delta_temperature_ignition, zrm_handle);
  if(status_other != ZERORK_STATUS_SUCCESS) error_state += 1;
#ifdef USE_OMP
  }
#endif
  if(status_mech != ZERORK_STATUS_SUCCESS) {
    printf("ERROR: Failed to parse mechanism file.\n");
    fflush(stdout); exit(-1);
  }
  if(status_options != ZERORK_STATUS_SUCCESS) {
    printf("WARNING: Failed to parse options file, continuing with default options.\n");
    fflush(stdout);
  }
  if(error_state>0) {
    printf("WARNING: Failed to set some cfd-plugin option(s).\n");
    fflush(stdout);
  }

  const char* cklogfilename = zerork::utilities::null_filename; //We already parsed in reactor manager no need to have another log
  // TODO: Avoid parsing twice/having two mechanisms
  zerork::mechanism mech(mechfilename, thermfilename, cklogfilename);


  int nSpc=mech.getNumSpecies();
  int nState=nSpc+1;
  int nStep=mech.getNumSteps();
  int nSpcStride = nSpc;

  int nReactors = inputFileDB.n_reactors();

  //Set up reactor initial states
  int nReactorsAlloc = nReactors;
  std::vector<double> reactorT(nReactors);
  std::vector<double> reactorP(nReactors);
  std::vector<double> reactorMassFrac(nReactors*nSpc);
  std::vector<double> reactorDPDT(nReactors,inputFileDB.dpdt());
  std::vector<double> reactorCost(nReactors);
  std::vector<double> reactorGpu(nReactors);
  std::vector<double> reactorESRC(nReactors, inputFileDB.e_src());
  std::vector<double> reactorYsrc(nReactors*nSpc, inputFileDB.y_src());
  std::vector<double> reactorIDT(nReactors,-1.0);
  std::vector<int> log_species_indexes(0);
  std::vector<std::string> log_species_names = inputFileDB.log_species();
  int n_print_reactors = std::min(inputFileDB.n_print_reactors(),nReactors);

  if( inputFileDB.log_species().size() > 0 ) {
     for (int k = 0; k < inputFileDB.log_species().size(); ++k) {
        std::string sp = inputFileDB.log_species()[k];
        int idx=mech.getIdxFromName(sp.c_str());
        if(idx == -1) {
          printf("WARNING: log_species %s not found in mechanism. Ignoring\n",sp.c_str());
        } else {
          log_species_indexes.push_back(idx);
        }
     }
  }
  if(inputFileDB.state_files().size() != 0) {
    //Load initial state from files
    if(inputFileDB.state_files().size() != nReactors) {
      printf("stateFiles present but not matching nReactors. Quitting.\n");
      exit(-1);
    }
    for(int k=0; k<nReactors; ++k) {
      FILE* stateFile = fopen(inputFileDB.state_files()[k].c_str(),"r");
      if(stateFile == NULL) {
        printf("Unable to read file: %s\n",inputFileDB.state_files()[k].c_str());
        exit(-1);
      }
      double val;
      for(int j=0;j<nSpc;++j) {
        int matches = fscanf(stateFile,"%lf",&val);
        if( matches != 1 ) {
          printf("Failed to find value in file %s\n",
                 inputFileDB.state_files()[k].c_str());
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
  } else if(inputFileDB.state_files_cfd().size() != 0) {
    int kReactor = 0;
    //Load initial state from files
    for(int k=0; k<inputFileDB.state_files_cfd().size() ; ++k) {
      FILE* stateFile = fopen(inputFileDB.state_files_cfd()[k].c_str(),"r");
      if(stateFile == NULL) {
        printf("Unable to read file: %s\n",inputFileDB.state_files_cfd()[k].c_str());
        exit(-1);
      }
      while(true) {
        if(kReactor >= nReactorsAlloc) {
          nReactorsAlloc += 100;
          reactorT.resize(nReactorsAlloc);
          reactorP.resize(nReactorsAlloc);
          reactorMassFrac.resize(nReactorsAlloc*nSpc);
          reactorDPDT.resize(nReactorsAlloc,inputFileDB.dpdt());
          reactorCost.resize(nReactorsAlloc,0);
          reactorGpu.resize(nReactorsAlloc,0);
          reactorIDT.resize(nReactorsAlloc,-1.0);
          reactorESRC.resize(nReactorsAlloc, inputFileDB.e_src());
          reactorYsrc.resize(nReactorsAlloc*nSpc, inputFileDB.y_src());
        }
        int ival;
        double dval;
        // read index
        // N.B. index is ignored
        int matches = fscanf(stateFile,"%d",&ival);
        if( matches != 1 ) {
          break;
        }
        fscanf(stateFile,"%lf",&dval);
        reactorT[kReactor] = dval;
        fscanf(stateFile,"%lf",&dval);
        reactorP[kReactor] = dval;
        fscanf(stateFile,"%d",&ival);
        reactorCost[kReactor] = (double) ival;
        fscanf(stateFile,"%lf",&dval);
        reactorGpu[kReactor] = dval;
        for(int j=0;j<nSpc;++j) {
          matches = fscanf(stateFile,"%lf",&dval);
          if( matches != 1 ) {
            printf("Failed to find value in file %s\n",
                   inputFileDB.state_files()[k].c_str());
            exit(-1);
          }
          reactorMassFrac[kReactor*nSpc+j] = dval;
        }
        kReactor += 1;
      }
      fclose(stateFile);
    }
    nReactors = kReactor;
    reactorT.resize(nReactors);
    reactorP.resize(nReactors);
    reactorMassFrac.resize(nReactors*nSpc);
    reactorDPDT.resize(nReactors);
    reactorCost.resize(nReactors);
    reactorGpu.resize(nReactors);
    reactorIDT.resize(nReactors);
    reactorESRC.resize(nReactors);
    reactorYsrc.resize(nReactors*nSpc);
    printf("Read %d reactors from %d files\n",nReactors, inputFileDB.state_files_cfd().size());
  } else {
    //Set initial states using linear ranges of T,P,phi, and egr
    double phiMin, phiMax, egrMin, egrMax, TMin, TMax, pMin, pMax;
    phiMin = inputFileDB.phi_min();
    phiMax = inputFileDB.phi_max();
    egrMin = inputFileDB.egr_min();
    egrMax = inputFileDB.egr_max();
    TMin = inputFileDB.temperature_min();
    TMax = inputFileDB.temperature_max();
    pMin = inputFileDB.pressure_min();
    pMax = inputFileDB.pressure_max();

    reactorT = linspace(TMin,TMax,nReactors);
    reactorP = linspace(pMin,pMax,nReactors);
    std::vector<double> reactorPhi = linspace(phiMin,phiMax,nReactors);
    std::vector<double> reactorEgr = linspace(egrMin,egrMax,nReactors);

    std::vector<double> initFuelMassFrac(nSpc);
    std::vector<double> initOxidMassFrac(nSpc);
    getFracsFromCompMap(mech,inputFileDB.fuel_composition(),&initFuelMassFrac[0]);
    getFracsFromCompMap(mech,inputFileDB.oxidizer_composition(),&initOxidMassFrac[0]);

    setReactorMassFrac(nReactors, mech, &reactorMassFrac[0],
                       &initFuelMassFrac[0], &initOxidMassFrac[0],
                       &reactorPhi[0], &reactorEgr[0]);
    if(log_species_indexes.size() == 0) {
        std::map<std::string, double>::const_iterator mapit;
        std::map<std::string, double> comp = inputFileDB.fuel_composition();
        for(mapit = comp.begin(); mapit != comp.end(); ++mapit)
        {
            std::string spcName = mapit->first;
            int spcIdx = mech.getIdxFromName(spcName.c_str());
            //N.B. We already checked fuel names in getFracsFromCompMap;
            assert(spcIdx != -1);
            log_species_indexes.push_back(spcIdx);
            log_species_names.push_back(spcName);
        }
        comp = inputFileDB.oxidizer_composition();
        for(mapit = comp.begin(); mapit != comp.end(); ++mapit)
        {
            std::string spcName = mapit->first;
            int spcIdx = mech.getIdxFromName(spcName.c_str());
            //N.B. We already checked oxid names in getFracsFromCompMap;
            assert(spcIdx != -1);
            log_species_indexes.push_back(spcIdx);
            log_species_names.push_back(spcName);
        }
    }
  }

  double tend = inputFileDB.solution_time();
  int n_steps = inputFileDB.n_steps();

  // timing data
  double startTime,stopTime,otherTime,simTime;
  startTime=getHighResolutionTime();

  int rank = 0;
  int nranks = 1;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nranks);
#endif
  std::vector<std::shared_ptr<std::ofstream>> reactor_log_files(0);
  if(rank != 0) {
    nReactors = 0;
  } else {
    for( int k = 0; k < n_print_reactors; ++k ) {
       std::stringstream filename;
       filename << inputFileDB.reactor_history_file_prefix();
       filename << "_" << std::setfill('0') << std::setw(3) << k << ".hist";
       reactor_log_files.push_back(std::make_shared<std::ofstream>(filename.str()));
    }
  }

  if(rank == 0 && nranks > 1 && !inputFileDB.batched()) {
    printf("WARNING: nranks > 1 can not be used with option \"batched\" disabled.\n");
    printf("         Running in batched mode.\n");
  }

  zerork_status_t flag = ZERORK_STATUS_SUCCESS;
  int num_solution_failures = 0;
  double t = 0;
  double dt= tend/n_steps;
  for(int i = 0; i < n_steps; ++i) {
      if(!inputFileDB.batched() && nranks == 1) {
#ifdef USE_OMP
        #pragma omp parallel for reduction(+:num_solution_failures)
#endif
        for(int k = 0; k < nReactors; ++k) {
            ud.nsteps = 0;
            ud.time = 0.0;
            if(delta_temperature_ignition > 0) {
                zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_IGNITION_TIME, &reactorIDT[k], zrm_handle);
            }
            if(inputFileDB.app_owns_aux_fields()) {
                zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_COST, &reactorCost[k], zrm_handle);
                zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_GPU, &reactorGpu[k], zrm_handle);
            }
            if(inputFileDB.dpdt() != 0.0) {
                zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_DPDT, &reactorDPDT[k], zrm_handle);
            }
            if(inputFileDB.e_src() != 0.0) {
                zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_E_SRC, &reactorESRC[k], zrm_handle);
            }
            if(inputFileDB.y_src() != 0.0) {
                zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_Y_SRC, &reactorYsrc[k*nSpc], zrm_handle);
            }
            flag = zerork_reactor_solve(i, t, dt, 1, &reactorT[k], &reactorP[k],
                                        &reactorMassFrac[k*nSpc], zrm_handle);
            if(flag != ZERORK_STATUS_SUCCESS) num_solution_failures+=1;
        }
      } else {
        if(delta_temperature_ignition > 0) {
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_IGNITION_TIME, &reactorIDT[0], zrm_handle);
        }
        if(inputFileDB.app_owns_aux_fields()) {
            //For AMR/moving mesh codes
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_COST, &reactorCost[0], zrm_handle);
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_GPU, &reactorGpu[0], zrm_handle);
        }
        if(inputFileDB.dpdt() != 0) {
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_DPDT, &reactorDPDT[0], zrm_handle);
        }
        if(inputFileDB.e_src() != 0.0) {
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_E_SRC, &reactorESRC[0], zrm_handle);
        }
        if(inputFileDB.y_src() != 0.0) {
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_Y_SRC, &reactorYsrc[0], zrm_handle);
        }
        flag = zerork_reactor_solve(i, t, dt, nReactors, &reactorT[0], &reactorP[0],
                                    &reactorMassFrac[0], zrm_handle);
        if(flag != ZERORK_STATUS_SUCCESS) num_solution_failures+=1;
      }
      if(rank==0) {
          log_output(i, t+dt, n_print_reactors, reactorT, reactorP, reactorDPDT,
                     reactorMassFrac, reactorCost, reactorGpu,
                     log_species_indexes, log_species_names,
                     nSpc, reactor_log_files);
      }
      if(num_solution_failures != 0) {
        printf("Zero-RK CFD Plugin failed to solve.\n");
        break;
      }
      t += dt;
  }
  if(rank == 0 && delta_temperature_ignition > 0) {
    for(int k = 0; k < nReactors; ++k) {
      if(reactorIDT[k] > 0) {
        printf("reactor[%d] IDT: %g s\n", k, reactorIDT[k]);
      }
    }
  }

  stopTime=getHighResolutionTime();
  simTime=stopTime-startTime;
  printf("simTime : %g s\n",simTime);
#ifdef USE_OMP
#pragma omp parallel
  {
#endif
  zerork_reactor_free(zrm_handle);
#ifdef USE_OMP
  }
#endif
}

int main(int argc, char **argv)
{
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
  zerork_reactor(argc, argv);
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
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

static void log_output(int step, double time, int n_print_reactors,
                       std::vector<double> reactorT, std::vector<double> reactorP,
                       std::vector<double> reactorDPDT, std::vector<double> reactorMassFrac,
                       std::vector<double> reactorCost, std::vector<double> reactorGpu,
                       std::vector<int> log_species_indexes, std::vector<std::string> log_species_names,
                       int nsp, std::vector<std::shared_ptr<std::ofstream>> reactor_log_files)
{
      int n_log_species = log_species_indexes.size();
      for (int k = 0; k < n_print_reactors; ++k) {
            std::ofstream& rlf = *reactor_log_files[k];
            if(step == 0) {
              rlf << "#";
              rlf << std::setw(12) <<  "step";
              rlf << std::setw(17) <<  "time";
              rlf << std::setw(17) <<  "temperature";
              rlf << std::setw(17) <<  "pressure";

              for(int j = 0; j < n_log_species; ++j) {
                  rlf << std::setw(17) << log_species_names[j];
              }
              rlf << std::endl;
            }

            rlf << std::setw(13) <<  step;
            rlf << std::setw(17) <<  time;
            rlf << std::setw(17) <<  reactorT[k];
            rlf << std::setw(17) <<  reactorP[k];

            for(int j = 0; j < n_log_species; ++j) {
                int spcIdx = log_species_indexes[j];
                double mf = reactorMassFrac[k*nsp+spcIdx];
                rlf << std::setw(17) << mf;
            }
            rlf << std::endl;
      }
}


