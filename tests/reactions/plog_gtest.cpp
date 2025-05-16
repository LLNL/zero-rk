#include <math.h>
#include <vector>
#include <cstdlib>

#include <zerork/mechanism.h>

#include <gtest/gtest.h>

// ---------------------------------------------------------------------------
// test constants
// ---------------------------------------------------------------------------
static const double OK_DOUBLE = 6.0e-15; // acceptable relative tolerance

static const char MECH_FILENAME[]  = "mechanisms/ideal/plog_test.mech";
static const char THERM_FILENAME[] = "mechanisms/ideal/const_specific_heat.therm";
static const char PARSER_LOGNAME[] = "parser.log";

static const int NUM_SPECIES = 2;
static const int NUM_REACTIONS = 8;

static const double LOCAL_GAS_CONST = 8.314462618e3;
static const double ATOMIC_A0 =  2.5;
static const double ATOMIC_A5 =  2.54737667e+04;
static const double ATOMIC_A6 = -4.46703359e-01; 
static const double DIATOMIC_A0 =  3.5;
static const double DIATOMIC_A5 = -1.04352500e+03;
static const double DIATOMIC_A6 = -4.22439182e+00;

static double GetEquilibriumConst(const double temperature,
                                  const double pressure);

static double GetRateConstant(zerork::mechanism *mech,
                              const double temperature,
                              const double pressure,
                              const std::vector<double> &mole_fractions,
                              const int reaction_id,
                              const int reaction_dir);
static void GetUniformComposition(zerork::mechanism *mech,
                                  std::vector<double> *composition);

static double kfwd_reaction07(const double p, const double T);
static const double KP1atm=1.01325e5;
static const double KEnergyToTemp=4184.0/LOCAL_GAS_CONST;  // [K/(cal/mol)]
static const double KConcentrationToKMol = 1000.0; // [mol/cm^3 -> kmol/m^3]

static bool NearScalar(const double a,
                      const double b);
static bool NearVector(size_t n,
                       const double a[],
                       const double b[]);

// ---------------------------------------------------------------------------
// test fixture for the Zero-RK mechanism 
class MechanismTestFixture: public ::testing::Test 
{ 
 public: 
  MechanismTestFixture( ) { 
    // initialization code here
    mechanism_    = NULL;
    const char * ZERORK_DATA_DIR = std::getenv("ZERORK_DATA_DIR");
    std::string load_mech(MECH_FILENAME);
    std::string load_therm(THERM_FILENAME);
    if(ZERORK_DATA_DIR == nullptr) {
      load_mech = std::string("../../") + load_mech;
      load_therm = std::string("../../") + load_therm;
    } else {
      load_mech = std::string(ZERORK_DATA_DIR) + "/" + load_mech;
      load_therm = std::string(ZERORK_DATA_DIR) + "/" + load_therm;
    }
    mechanism_ = new zerork::mechanism(load_mech.c_str(),
                                       load_therm.c_str(),
                                       PARSER_LOGNAME);
  } 

  void SetUp( ) { 
    // code here will execute just before the test ensues 
  }

  void TearDown( ) { 
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be
  }

  ~MechanismTestFixture( )  { 
    // cleanup any pending stuff, but no exceptions allowed
    if(mechanism_ != NULL) {
      delete mechanism_;
    }
  }

  // put in any custom data members that you need 
  zerork::mechanism    *mechanism_;
};

TEST_F (MechanismTestFixture, Allocation)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";
  
  const int num_species   = mechanism_->getNumSpecies();
  const int num_reactions = mechanism_->getNumReactions();

  ASSERT_TRUE(NUM_SPECIES == num_species) <<
    "mechanism_->getNumSpecies()";

  ASSERT_TRUE(NUM_REACTIONS == num_reactions) <<
    "mechanism_->getNumReactions()";
}

TEST_F (MechanismTestFixture, Reaction01)
{

  ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";
  
  const int num_species   = mechanism_->getNumSpecies();
  const int num_reactions = mechanism_->getNumReactions();
  std::vector<double> mole_fractions;
  double K_forward, K_reverse;
  double temperature, pressure;
  
  ASSERT_TRUE(NUM_SPECIES == num_species) <<
    "mechanism_->getNumSpecies()";

  ASSERT_TRUE(NUM_REACTIONS == num_reactions) <<
    "mechanism_->getNumReactions()";

  GetUniformComposition(mechanism_,&mole_fractions);

  temperature = 1000.0;
  pressure = 1.01325e5;
  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              0,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)
  K_reverse = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              0,   // reaction_id
                              -1); // reaction direction (fwd >= 0, rev < 0)
                      
  EXPECT_TRUE(NearScalar(K_forward, 
                         1.0e6*pow(temperature,2.0)*exp(-1000.0/temperature))) << "K_forward - single pressure plog at  1 atm (exact table point)";
  EXPECT_TRUE(NearScalar(K_reverse, 0.0)) << "K_reverse - single pressure plog at pressure point";

  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              0.1*pressure,
                              mole_fractions,
                              0,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)
                      
  EXPECT_TRUE(NearScalar(K_forward, 
                         1.0e6*pow(temperature,2.0)*exp(-1000.0/temperature))) << "K_forward - single pressure plog at 0.1 atm (under table)";

  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              10.0*pressure,
                              mole_fractions,
                              0,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)
                      
  EXPECT_TRUE(NearScalar(K_forward, 
                         1.0e6*pow(temperature,2.0)*exp(-1000.0/temperature))) << "K_forward - single pressure plog at 10 atm (over table)";

}


TEST_F (MechanismTestFixture, Reaction02)
{
    ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";
  
  const int num_species   = mechanism_->getNumSpecies();
  const int num_reactions = mechanism_->getNumReactions();
  std::vector<double> mole_fractions;
  double K_forward, K_reverse;
  double temperature, pressure;
  
  ASSERT_TRUE(NUM_SPECIES == num_species) <<
    "mechanism_->getNumSpecies()";

  ASSERT_TRUE(NUM_REACTIONS == num_reactions) <<
    "mechanism_->getNumReactions()";

  GetUniformComposition(mechanism_,&mole_fractions);

  temperature = 800.0;
  pressure = 1.01325e6;

  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              2.0*pressure,
                              mole_fractions,
                              1,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)
                      
  EXPECT_TRUE(NearScalar(K_forward, 
                         -1.0e6*pow(temperature,2.0)*exp(-1000.0/temperature))) << "K_forward - single pressure plog with negative A-Factor";

}

TEST_F (MechanismTestFixture, Reaction03)
{
    ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";
  
  const int num_species   = mechanism_->getNumSpecies();
  const int num_reactions = mechanism_->getNumReactions();
  std::vector<double> mole_fractions;
  double K_forward, K_reverse;
  double temperature, pressure;
  
  ASSERT_TRUE(NUM_SPECIES == num_species) <<
    "mechanism_->getNumSpecies()";

  ASSERT_TRUE(NUM_REACTIONS == num_reactions) <<
    "mechanism_->getNumReactions()";

  GetUniformComposition(mechanism_,&mole_fractions);

  temperature = 1500;
  pressure = 1.0e5;

  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              2.0*pressure,
                              mole_fractions,
                              2,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)
                      
  EXPECT_TRUE(NearScalar(K_forward, 
                         9.0e5*pow(temperature,2.0)*exp(-1000.0/temperature))) << "K_forward - single pressure plog with multiple Arrhenius lines";

}

TEST_F (MechanismTestFixture, Reaction04)
{
    ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";
  
  const int num_species   = mechanism_->getNumSpecies();
  const int num_reactions = mechanism_->getNumReactions();
  std::vector<double> mole_fractions;
  double K_forward, K_reverse;
  double temperature, pressure;
  
  ASSERT_TRUE(NUM_SPECIES == num_species) <<
    "mechanism_->getNumSpecies()";

  ASSERT_TRUE(NUM_REACTIONS == num_reactions) <<
    "mechanism_->getNumReactions()";

  GetUniformComposition(mechanism_,&mole_fractions);

  temperature = 3200; // high temprature avoids very large equilibrium
  pressure = 1.0e6;   // values that create a reverse rate constant
                      // that is less than 

  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              3,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)

  K_reverse = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              3,  // reaction_id
                              -1); // reaction direction (fwd >= 0, rev < 0)

  double K_equilibrium = K_forward/K_reverse;
  double K_equilibrium_check = GetEquilibriumConst(temperature,
                                                   pressure);
  double K_forward_check = 1.0e6*pow(temperature,2.0)*exp(-1000.0/temperature);
  double K_reverse_check = K_forward_check/K_equilibrium_check;

  EXPECT_TRUE(NearScalar(K_forward, K_forward_check)) << "K_forward - single pressure plog";
  //printf("Keq = %.18g\n",K_equilibrium);

  EXPECT_TRUE(NearScalar(K_reverse, K_reverse_check)) << "K_reverse - calculated from equilibrium const for a single pressure plog";

}

TEST_F (MechanismTestFixture, Reaction05)
{
    ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";
  
  const int num_species   = mechanism_->getNumSpecies();
  const int num_reactions = mechanism_->getNumReactions();
  std::vector<double> mole_fractions;
  double K_forward, K_reverse;
  double temperature, pressure;
  
  ASSERT_TRUE(NUM_SPECIES == num_species) <<
    "mechanism_->getNumSpecies()";

  ASSERT_TRUE(NUM_REACTIONS == num_reactions) <<
    "mechanism_->getNumReactions()";

  GetUniformComposition(mechanism_,&mole_fractions);

  temperature = 3600; // high temprature avoids very large equilibrium
  pressure = 1.0e4;   // values that create a reverse rate constant
                      // that is less than 

  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              4,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)

  K_reverse = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              4,  // reaction_id
                              -1); // reaction direction (fwd >= 0, rev < 0)

  double K_equilibrium = K_forward/K_reverse;
  double K_equilibrium_check = GetEquilibriumConst(temperature,
                                                   pressure);
  double K_forward_check = -1.0e6*pow(temperature,2.0)*exp(-1000.0/temperature);
  double K_reverse_check = K_forward_check/K_equilibrium_check;

  EXPECT_TRUE(NearScalar(K_forward, K_forward_check)) << "K_forward - single pressure plog";
  //printf("Keq = %.18g\n",K_equilibrium);

  EXPECT_TRUE(NearScalar(K_reverse, K_reverse_check)) << "K_reverse - calculated from equilibrium const for a single pressure plog";

}

TEST_F (MechanismTestFixture, Reaction06)
{
    ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";
  
  const int num_species   = mechanism_->getNumSpecies();
  const int num_reactions = mechanism_->getNumReactions();
  std::vector<double> mole_fractions;
  double K_forward, K_reverse;
  double temperature, pressure, power_law;
  
  ASSERT_TRUE(NUM_SPECIES == num_species) <<
    "mechanism_->getNumSpecies()";

  ASSERT_TRUE(NUM_REACTIONS == num_reactions) <<
    "mechanism_->getNumReactions()";

  GetUniformComposition(mechanism_,&mole_fractions);

  temperature = 3600;
  pressure = 1.0e5;
  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              5,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)

  EXPECT_TRUE(NearScalar(K_forward, 6.2)) << "K_forward - p = 1 bar (3 table pressures with two Arrhenius lines each)";

  pressure = 1.01325e5;
  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              5,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)

  EXPECT_TRUE(NearScalar(K_forward, 6.2)) << "K_forward - p = 1 atm (3 table pressures with two Arrhenius lines each)";

  pressure = 1.5*1.01325e5;
  power_law = log(1.3/6.2)/log(2.0/1.0);
  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              5,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)

  EXPECT_TRUE(NearScalar(K_forward, 6.2*pow(1.5, power_law))) << "K_forward - p = 1.5 atm (3 table pressures with two Arrhenius lines each)";

  pressure = 2.0*1.01325e5;
  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              5,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)

  EXPECT_TRUE(NearScalar(K_forward, 1.3)) << "K_forward - p = 2 atm (3 table pressures with two Arrhenius lines each)";

  pressure = 3.5*1.01325e5;
  power_law = log(10.01/1.3)/log(4.0/2.0);
  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              5,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)

  EXPECT_TRUE(NearScalar(K_forward, 1.3*pow(3.5/2.0, power_law))) << "K_forward - p = 3.5 atm (3 table pressures with two Arrhenius lines each)";

  pressure = 4.0*1.01325e5;
  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              5,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)

  EXPECT_TRUE(NearScalar(K_forward,10.01)) << "K_forward - p = 4 atm (3 table pressures with two Arrhenius lines each)";

  pressure = 5.0*1.01325e5;
  K_forward = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              5,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)

  EXPECT_TRUE(NearScalar(K_forward,10.01)) << "K_forward - p = 5 atm (3 table pressures with two Arrhenius lines each)";

}

TEST_F (MechanismTestFixture, Reaction07)
{
    ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";
  
  const int num_species   = mechanism_->getNumSpecies();
  const int num_reactions = mechanism_->getNumReactions();
  std::vector<double> mole_fractions;
  double K_forward, K_forward_check;
  double temperature, pressure, pressure_multiplier;
  
  ASSERT_TRUE(NUM_SPECIES == num_species) <<
    "mechanism_->getNumSpecies()";

  ASSERT_TRUE(NUM_REACTIONS == num_reactions) <<
    "mechanism_->getNumReactions()";

  GetUniformComposition(mechanism_,&mole_fractions);

  pressure = 1.0e3;
  temperature = 500.0;

  while(pressure < 2.0e7) {

    K_forward = GetRateConstant(mechanism_,
                                temperature,
                                pressure,
                                mole_fractions,
                                6,  // reaction_id
                                1); // reaction direction (fwd >= 0, rev < 0)
    K_forward_check = kfwd_reaction07(pressure, temperature);

    EXPECT_TRUE(NearScalar(K_forward,K_forward_check)) << "K_forward at p = " 
					               << pressure << ", T = "
                                                       << temperature;
    pressure *= 2.0;
    temperature += 200.0;
  }
}

TEST_F (MechanismTestFixture, Reaction08)
{
    ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";
  
  const int num_species   = mechanism_->getNumSpecies();
  const int num_reactions = mechanism_->getNumReactions();
  std::vector<double> mole_fractions;
  double K_forward, K_forward_check;
  double temperature, pressure, pressure_multiplier;
  
  ASSERT_TRUE(NUM_SPECIES == num_species) <<
    "mechanism_->getNumSpecies()";

  ASSERT_TRUE(NUM_REACTIONS == num_reactions) <<
    "mechanism_->getNumReactions()";

  GetUniformComposition(mechanism_,&mole_fractions);

  pressure = 1.0e3;
  temperature = 500.0;

  while(pressure < 2.0e7) {

    K_forward = GetRateConstant(mechanism_,
                                temperature,
                                pressure,
                                mole_fractions,
                                7,  // reaction_id
                                1); // reaction direction (fwd >= 0, rev < 0)
    K_forward_check = kfwd_reaction07(pressure, temperature);

    EXPECT_TRUE(NearScalar(K_forward,K_forward_check)) << "K_forward at p = " 
					               << pressure << ", T = "
                                                       << temperature;
    pressure *= 2.0;
    temperature += 200.0;
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


static double GetRateConstant(zerork::mechanism *mech,
                              const double temperature,
                              const double pressure,
                              const std::vector<double> &mole_fractions,
                              const int reaction_id,
                              const int reaction_dir)
{
  const double gas_constant = mech->getGasConstant();
  const int num_species = mech->getNumSpecies();
  const int num_reactions = mech->getNumReactions();
  double Krxn;

  std::vector<double> concentrations;
  std::vector<double> forward_rate_coefficients;
  std::vector<double> reverse_rate_coefficients;

  concentrations.clear();
  concentrations.assign(num_species, 0.0);
  forward_rate_coefficients.clear();
  forward_rate_coefficients.assign(num_reactions, 0.0);
  reverse_rate_coefficients.clear();
  reverse_rate_coefficients.assign(num_reactions, 0.0);

  for(int j=0; j<num_species; ++j) {
    concentrations[j] = mole_fractions[j]*pressure/(gas_constant*temperature);
  }
  mech->getKrxnFromTC(temperature,
                      &concentrations[0],
                      &forward_rate_coefficients[0],
                      &reverse_rate_coefficients[0]);
  Krxn = ((reaction_dir < 0) ? reverse_rate_coefficients[reaction_id] :
          forward_rate_coefficients[reaction_id]);

  return Krxn;

}
static void GetUniformComposition(zerork::mechanism *mech,
                                  std::vector<double> *composition)
{
  const int num_species = mech->getNumSpecies();
  composition->clear();
  composition->assign(num_species,1.0/(double)num_species);
}

static double GetEquilibriumConst(const double temperature,
                                  const double pressure)
{
  // A + A <=> A2 (diatomic)
  // S/Ru
  double   atomic_entropy =   ATOMIC_A0*log(temperature) +   ATOMIC_A6;
  double diatomic_entropy = DIATOMIC_A0*log(temperature) + DIATOMIC_A6;
  // H/RuT
  double   atomic_enthalpy =   ATOMIC_A0 +   ATOMIC_A5/temperature;
  double diatomic_enthalpy = DIATOMIC_A0 + DIATOMIC_A5/temperature;

  double delta_entropy  = -2.0*atomic_entropy  + diatomic_entropy;
  double delta_enthalpy = -2.0*atomic_enthalpy + diatomic_enthalpy;

  double Keq = pow(1.01325e5/(LOCAL_GAS_CONST*temperature), -1.0)*
    exp(delta_entropy-delta_enthalpy);

  return Keq;
}





static bool NearScalar(const double a,
                       const double b)
{
  double weight = 0.5*(fabs(a)+fabs(b));

  if(weight >= OK_DOUBLE*OK_DOUBLE) {
    // check the normalized difference
    if(fabs(a-b)/weight > OK_DOUBLE) {

      printf("#         %24.18e != %24.18e\n",a,b);
      return false;
    } 
  }
  return true;
}

static bool NearVector(size_t n,
                       const double a[],
                       const double b[])
{
  bool near_scalar;
  for(size_t j=0; j<n; ++j) {
    near_scalar = NearScalar(a[j],b[j]);
    if(near_scalar == false) {
      return false;
    }
  }
  return true;
}

// actual plog in use from tpgme mechanism
static double kfwd_reaction07(const double p, const double T)
{
  int id;
  double p_atm = p/KP1atm; // convert pressure to atmospheres

  const int n_pts = 5;
  double p_pts[]={0.01, 0.1, 1.0, 10.0, 100.0};
  double kfwd[n_pts];  
  double log_kfwd;

  // plog/ 0.0100   3.502e+005        1.441        -3244.0/
  kfwd[0] =         3.502e+005*pow(T, 1.441)*exp(-(-3244.0*KEnergyToTemp)/T);

  // plog/ 0.1000   8.854e+005        1.327        -2975.0/
  kfwd[1] =         8.854e+005*pow(T, 1.327)*exp(-(-2975.0*KEnergyToTemp)/T);

  // plog/ 1.0000   1.650e+007        0.973        -2010.0/
  kfwd[2] =         1.650e+007*pow(T, 0.973)*exp(-(-2010.0*KEnergyToTemp)/T);

  // plog/ 10.000   5.374e+009        0.287          280.0/
  kfwd[3] =         5.374e+009*pow(T, 0.287)*exp(-(  280.0*KEnergyToTemp)/T);

  // plog/ 100.00   9.494e+018       -2.199         9769.0/
  kfwd[4] =         9.494e+018*pow(T,-2.199)*exp(-( 9769.0*KEnergyToTemp)/T);  

  if(p_atm < p_pts[0]) {
    return kfwd[0]/KConcentrationToKMol; // no extrapolation
  }
  else if(p_atm >= p_pts[n_pts-1]) {
    return kfwd[4]/KConcentrationToKMol; // no extrapolation
  }
  else {
    id = -1;
    for(int j=0; j<(n_pts-1); ++j) {
      if(p_pts[j] <= p_atm && p_atm < p_pts[j+1]) {
        id = j;
        break;
      }
    }
    if(id == -1) {
      printf("ERROR: In kfwd_reaction05(...), failed to find pressure range\n");
      printf("       for p = %.18g [Pa] = %.18g [atm]\n",
             p,p_atm);
      exit(-1);
    }
  }
  log_kfwd = log(kfwd[id]) + 
    log(kfwd[id+1]/kfwd[id])*log(p_atm/p_pts[id])/log(p_pts[id+1]/p_pts[id]);
  // reaction is second order in the forward direction so the rate constant
  // must be converted from [cm^3/mol/s] to [m^3/kmol/s]
  return exp(log_kfwd)/KConcentrationToKMol;
}
