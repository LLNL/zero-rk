#include "soot.h"

// Defined constants from FlameMaster
static const double RGAS = 8314.462;    // [J / kmole K]
static const double AVOGADRO = 6.0221367e26;   // [1 / kmole]
static const double fMolarMassSoot = 12.01115; // [kg/kmole]
static const double fSootDensity   = 1800.0;   // [kg/m^3]
static const double Df             = 1.8;
static const double pi             = 3.1415926535897932385;
static const double CarbonToDiam   = pow(6.0*fMolarMassSoot/
                                         (pi*fSootDensity*AVOGADRO),
                                         1.0/3.0);
static const double CRed           = sqrt(0.5 * pi * RGAS / fMolarMassSoot)
                                     * CarbonToDiam*CarbonToDiam
                                     * AVOGADRO;

static const double nucl_nbrC = 40.0;
static const double nucl_surf = 11.69607095; // pow(nucl_nbrC,2.0/3.0);
static const double nucl_nbrH = 32.0;

// List of Species for Nucleation
int nArom;
int *AromInd;
double *AromNbrC;
double *AromNbrH;
double *AromStick;

// Local Dimer parameters
// -> to be updated at each grid points
double dimer_conc;
double dimer_rate;
double dimer_nbrC;
double dimer_nbrH;


//
void InitializePAH(void *user_data)
{
  FlameParams *params = (FlameParams *)user_data;

  double stick_multiplier = 1.5e-11;
  std::vector<string> state_name;

  int num_species = params->reactor_->GetNumSpecies();
  std::vector<int> num_hydrogen, num_carbon;
  num_hydrogen.assign(num_species, 0);
  num_carbon.assign(num_species, 0);

  params->reactor_->GetSpeciesHydrogenCount(&num_hydrogen[0]);
  params->reactor_->GetSpeciesCarbonCount(&num_carbon[0]);

  string prefix = "MassFraction_";
  string line;
  ifstream PAH_file;
  PAH_file.open("PAH_file");
  if(!PAH_file.is_open()) {
    cerr << "# Error opening PAH_file\n";
    exit(-1);
  }
  while(getline(PAH_file,line)) {
    line = prefix+line;
    state_name.push_back(line);
  }

  nArom = state_name.size();
  printf("# nArom = %d\n", nArom);

  AromInd   = new int[nArom];
  AromNbrC  = new double[nArom];
  AromNbrH  = new double[nArom];
  AromStick = new double[nArom]; //AromStic = 1.5e-11*MW^4 with MW in grams/mole

  for (int i=0; i<nArom; i++) {
    if ((AromInd[i] = params->reactor_->GetIdOfState(state_name[i].c_str()) )    != -1) {
      AromNbrC[i] = num_carbon[AromInd[i]];
      AromNbrH[i] = num_hydrogen[AromInd[i]];
      AromStick[i] = stick_multiplier*pow(AromNbrC[i]*12.0+AromNbrH[i], 4.0);
    }
    // Hard-code some values which differ from scaling above
    if(state_name[i]=="MassFraction_NAPH") {AromStick[i]=0.002;}
    if(state_name[i]=="MassFraction_A2R5") {AromStick[i]=0.004;}
    printf("# Aromatics species name: %s, ind: %d, nC: %g, nH: %g, stick: %g\n", state_name[i].c_str(), AromInd[i], AromNbrC[i], AromNbrH[i], AromStick[i]);
  }

}

double CollisionEfficiency(const int i,
                           const int j,
                           const double temperature,
                           const double density,
                           const double molecular_mass[])
{
  double nbrC_PAH = AromNbrC[i]+AromNbrC[j];
  double diam_PAH = CarbonToDiam * pow(nbrC_PAH,1.0/3.0);
  double mass_PAH = molecular_mass[AromInd[i]] + molecular_mass[AromInd[j]];

  static const double diam_Gas = 3.621e-10; // N2 sigma
  static const double mass_Gas = 28.0;      // N2 molecular weight
  static const double cEff = sqrt (0.5 * pi * RGAS) * AVOGADRO;

  double beta = cEff * sqrt(temperature*(1.0/mass_PAH+1.0/mass_Gas))
    * pow(diam_PAH+diam_Gas,2.0);

  double dt = 1.0 / (beta * density / mass_Gas);
  double tau  = (0.0435*mass_PAH - 7.0943) * 1e-12;
  double efficiency = 0.02*tau/dt;

  return efficiency;
}

void ComputeDimerProdRate(void *user_data,
                          const double state[],
                          double dimer_prod_rate[])
{
  FlameParams *params = (FlameParams *)user_data;
  const int num_species = params->reactor_->GetNumSpecies();
  const double temperature = state[num_species+1]*
    params->ref_temperature_;
  const double density = 1.0/state[num_species];

  std::vector<double> molecular_mass;
  molecular_mass.assign(num_species, 0.0);
  params->reactor_->GetSpeciesMolecularWeight(&molecular_mass[0]);

  // Reset
  dimer_prod_rate[0] = 0.0;
  dimer_prod_rate[1] = 0.0;
  dimer_prod_rate[2] = 0.0;

  for (int i=0; i<nArom; i++){
    double wDimer = 0.5 * GetBetaDimer(temperature,AromNbrC[i],AromNbrC[i]) *
      density * state[AromInd[i]] / molecular_mass[AromInd[i]] *
      density * state[AromInd[i]] / molecular_mass[AromInd[i]];
    if (wDimer<0.0) wDimer = 0.0;

    wDimer *= AromStick[i];
    dimer_prod_rate[0] += wDimer;               // number density
    dimer_prod_rate[1] += 2.0 * AromNbrC[i] * wDimer; // nbrC density
    dimer_prod_rate[2] += 2.0 * AromNbrH[i] * wDimer; // nbrH density
  }

  if (dimer_prod_rate[0]==0.0) dimer_prod_rate[0] = 1.0e-60;

}

double GetBetaDimer(const double temperature,
                    const double i,
                    const double j)
{
  double c = CRed*sqrt(temperature);
  double beta_dimer = c * sqrt(i+j) * pow(pow(i,1.0/3.0)+pow(j,1.0/3.0),2.0) / sqrt(i*j);

  return beta_dimer;
}

void CalcRhoDot(void *user_data,
                const double state[],
                double rho_dot)
{
  FlameParams *params = (FlameParams *)user_data;
  const int num_species = params->reactor_->GetNumSpecies();
  const double temperature = state[num_species+1]*
    params->ref_temperature_;
  const double density = 1.0/state[num_species];

  std::vector<double> molecular_mass;
  molecular_mass.assign(num_species, 0.0);
  params->reactor_->GetSpeciesMolecularWeight(&molecular_mass[0]);

  // Reset
  rho_dot = 0.0;

  for (int i=0; i<nArom; i++){
      double wDimer = 0.5 * GetBetaDimer(temperature,AromNbrC[i],AromNbrC[i]) *
        density * state[AromInd[i]] / molecular_mass[AromInd[i]] *
        density * state[AromInd[i]] / molecular_mass[AromInd[i]];
      if (wDimer<0.0) wDimer = 0.0;

      wDimer *= AromStick[i];

      // Rate of nucleation
      rho_dot -= 2.0 * molecular_mass[AromInd[i]] * wDimer;
  }

}

void UpdateProductionRates(void *user_data,
		           const double state[],
		           double prod_rate[])
{
  FlameParams *params = (FlameParams *)user_data;
  const int num_species = params->reactor_->GetNumSpecies();
  const double temperature = state[num_species+1]*
    params->ref_temperature_;
  const double density = 1.0/state[num_species];

  std::vector<double> molecular_mass;
  molecular_mass.assign(num_species, 0.0);
  params->reactor_->GetSpeciesMolecularWeight(&molecular_mass[0]);

  for (int i=0; i<nArom; i++){
    double wDimer = 0.5 * GetBetaDimer(temperature,AromNbrC[i],AromNbrC[i]) *
      density * state[AromInd[i]] / molecular_mass[AromInd[i]] *
      density * state[AromInd[i]] / molecular_mass[AromInd[i]];
    if (wDimer<0.0) wDimer = 0.0;

    wDimer *= AromStick[i];

    // Rate of nucleation
    prod_rate[AromInd[i]] -= 2.0*molecular_mass[AromInd[i]] * wDimer/density;
  }

}

void ComputeDimerParticles(void *user_data,
                           const double state[])
{
  FlameParams *params = (FlameParams *)user_data;
  const int num_species = params->reactor_->GetNumSpecies();
  const double temperature = state[num_species+1]*
    params->ref_temperature_;
    double dimer_prod_rate[3];

  ComputeDimerProdRate(params, &state[0], &dimer_prod_rate[0]);

  // Set number of carbon atoms
  dimer_rate = dimer_prod_rate[0];
  if (dimer_prod_rate[0]>1.0e-60){
    dimer_nbrC = dimer_prod_rate[1] / dimer_prod_rate[0];
    dimer_nbrH = dimer_prod_rate[2] / dimer_prod_rate[0];
  } else {
    dimer_nbrC = nucl_nbrC;
    dimer_nbrH = nucl_nbrH;
  }

  // Sink term due to "nucleation"
  double betaN = 0.0;
  betaN = GetBetaNucleation(temperature, dimer_nbrC);

  // Solve quadratic equation for fictive dimer species
  double Delta = 4.0*betaN*dimer_prod_rate[0];
  if (Delta>=0.0)
    dimer_conc = sqrt(Delta) / (2.0*betaN);
  else{
    cout << dimer_prod_rate[1] << "\t" << dimer_prod_rate[2] << "\n";
    cout << dimer_nbrC    << "\t" << dimer_nbrH    << "\n";
    cerr << "betaN: " << betaN << "\n";
    cerr << "dimer_prod_rate: " << dimer_prod_rate[0] << "\n";
    cerr << "error: negative Delta: " << Delta << "\n";
    exit(2);
  }

  // Check negative concentration
  if (dimer_conc < 0.0)
    dimer_conc = 0.0;

}

double GetBetaNucleation(const double temperature,
                         const double i)
{
    double c = CRed*sqrt(temperature);
    double beta_nucleation = 2.2*c*4.0*sqrt(2.0) * pow(i,1.0/6.0);
    return beta_nucleation;
}



void WriteDimerProdRate(void *user_data,
                        const double state[])
{
  FlameParams *params = (FlameParams *)user_data;
  const int num_states = params->reactor_->GetNumStates();
  const int num_local_points = params->num_local_points_;
  const int num_species = params->reactor_->GetNumSpecies();
  const int nover = params->nover_;
  int npes = params->npes_;
  int my_pe = params->my_pe_;
  MPI_Comm comm = params->comm_;
  double dimer_prod_rate[3];
  double local_sum = 0.0;
  double sum_dimer = 0.0;

  std::vector<double> dimer_prod_rate_local;
  dimer_prod_rate_local.assign(num_local_points, 0.0);
  std::vector<double> dimer_prod_rate_all;
  dimer_prod_rate_all.assign(num_local_points*npes, 0.0);

  for(int j=0; j<num_local_points; j++) {
    ComputeDimerProdRate(params,
                         &state[j*num_states],
                         &dimer_prod_rate[0]);
    dimer_prod_rate_local[j] = dimer_prod_rate[0];
  }

  long int dsize = num_local_points;
  int nodeDest = 0;

  MPI_Gather(&dimer_prod_rate_local[0],
             dsize,
             PVEC_REAL_MPI_TYPE,
             &dimer_prod_rate_all[0],
             dsize,
             PVEC_REAL_MPI_TYPE,
             nodeDest,
             comm);

  for(int j=0; j<num_local_points; ++j) {
    int jext = j + nover;
    local_sum += dimer_prod_rate_local[j]*state[j*num_states + num_species]/
      params->convection_velocity_[jext]*params->dzm_local_[jext];
  }
  MPI_Allreduce(&local_sum,&sum_dimer,1,PVEC_REAL_MPI_TYPE,MPI_SUM,comm);

  if(my_pe == 0) {
    FILE * dimerFile;
    dimerFile = fopen("dimer","w");
    fprintf(dimerFile,"# Integrated soot yield: %14.7e\n", sum_dimer);
    fprintf(dimerFile,"#Z      Dimer\n");
    for(int j=0; j<num_local_points*npes; j++) {
      fprintf(dimerFile,"%14.7e   %14.7e\n", params->z_[j], dimer_prod_rate_all[j]);
    }
  fclose(dimerFile);
  }

}


void UpdateDimerProdRate(void *user_data,
                        const double state[])
{
  FlameParams *params = (FlameParams *)user_data;
  const int num_states = params->reactor_->GetNumStates();
  const int num_local_points = params->num_local_points_;
  const int num_species = params->reactor_->GetNumSpecies();
  const int nover = params->nover_;
  int npes = params->npes_;
  MPI_Comm comm = params->comm_;
  double dimer_prod_rate[3];
  double local_sum = 0.0;
  double sum_dimer = 0.0;

  std::vector<double> dimer_prod_rate_local;
  dimer_prod_rate_local.assign(num_local_points, 0.0);
  std::vector<double> dimer_prod_rate_all;
  dimer_prod_rate_all.assign(num_local_points*npes, 0.0);

  for(int j=0; j<num_local_points; j++) {
    ComputeDimerProdRate(params,
                         &state[j*num_states],
                         &dimer_prod_rate[0]);
    dimer_prod_rate_local[j] = dimer_prod_rate[0];
  }

  long int dsize = num_local_points;
  int nodeDest = 0;

  MPI_Gather(&dimer_prod_rate_local[0],
             dsize,
             PVEC_REAL_MPI_TYPE,
             &dimer_prod_rate_all[0],
             dsize,
             PVEC_REAL_MPI_TYPE,
             nodeDest,
             comm);

  for(int j=0; j<num_local_points; ++j) {
    int jext = j + nover;
    local_sum += dimer_prod_rate_local[j]*state[j*num_states + num_species]/
      params->convection_velocity_[jext]*params->dzm_local_[jext];
  }
  MPI_Allreduce(&local_sum,&sum_dimer,1,PVEC_REAL_MPI_TYPE,MPI_SUM,comm);

  params->Y_sootmax = sum_dimer;

}
