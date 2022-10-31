#include <math.h>
#include <vector>
#include <cstdlib>

#include <zerork/mechanism.h>

#include <gtest/gtest.h>

// ---------------------------------------------------------------------------
// test constants
// ---------------------------------------------------------------------------
static const double OK_DOUBLE = 6.0e-15; // acceptable relative tolerance

static const char MECH_FILENAME[]  = "mechanisms/ideal/troe_test.mech";
static const char THERM_FILENAME[] = "mechanisms/ideal/const_specific_heat.therm";
static const char PARSER_LOGNAME[] = "parser.log";

static const int NUM_SPECIES = 2;
static const int NUM_REACTIONS = 2;

static const double LOCAL_GAS_CONST = 8314.4598;
static const double ATOMIC_A0 =  2.5;
static const double ATOMIC_A5 =  2.54737667e+04;
static const double ATOMIC_A6 = -4.46703359e-01; 
static const double DIATOMIC_A0 =  3.5;
static const double DIATOMIC_A5 = -1.04352500e+03;
static const double DIATOMIC_A6 = -4.22439182e+00;

static double GetEquilibriumConst(const double temperature,
                                  const double pressure);

static double GetConcentrationSum(zerork::mechanism *mech,
                              const double temperature,
                              const double pressure,
                              const std::vector<double> &mole_fractions);
static double GetRateConstant(zerork::mechanism *mech,
                              const double temperature,
                              const double pressure,
                              const std::vector<double> &mole_fractions,
                              const int reaction_id,
                              const int reaction_dir);
static void GetUniformComposition(zerork::mechanism *mech,
                                  std::vector<double> *composition);

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
  double Khigh, Klow, K_forward_check, Fcenter;
  double Pr, fTerm, nTerm, log_10_Pr, Pcorr, Csum;
  
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

  Csum = GetConcentrationSum(mechanism_,
                             temperature,
                             pressure,
                             mole_fractions);
  Khigh = 2.0e12;  
  Klow = 2.0e6; //N.B. convertC
  Pr = Klow*Csum/Khigh;
  log_10_Pr=log10(Pr);

  Fcenter=(1.0-0.68)*exp(-temperature/78.0)
              +0.68*exp(-temperature/1995.0);
  fTerm=log10(Fcenter);
  nTerm=0.75-1.27*fTerm;
  log_10_Pr-=(0.4+0.67*fTerm);                // log10(Pr) + c
  log_10_Pr=log_10_Pr/(nTerm-0.14*log_10_Pr); // d = 0.14
  log_10_Pr*=log_10_Pr;
  fTerm/=(1.0+log_10_Pr);
  fTerm=pow(10.0,fTerm);
  Pcorr = fTerm*Pr/(1.0+Pr);

  K_forward_check = Khigh*Pcorr;
  EXPECT_TRUE(NearScalar(K_forward,K_forward_check)) << "K_forward - three term troe";
  EXPECT_TRUE(NearScalar(K_reverse, 0.0)) << "K_reverse - three term troe";

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
  double Khigh, Klow, K_forward_check, Fcenter;
  double Pr, fTerm, nTerm, log_10_Pr, Pcorr, Csum;
  
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
                              1,  // reaction_id
                              1); // reaction direction (fwd >= 0, rev < 0)
  K_reverse = GetRateConstant(mechanism_,
                              temperature,
                              pressure,
                              mole_fractions,
                              1,   // reaction_id
                              -1); // reaction direction (fwd >= 0, rev < 0)

  Csum = GetConcentrationSum(mechanism_,
                             temperature,
                             pressure,
                             mole_fractions);
  Khigh = 2.0e12;  
  Klow = 2.0e6; //N.B. convertC
  Pr = Klow*Csum/Khigh;
  log_10_Pr=log10(Pr);

  Fcenter=(1.0-0.68)*exp(-temperature/78.0)
              +0.68*exp(-temperature/1995.0);
  Fcenter += exp(-5590/temperature);  
  fTerm=log10(Fcenter);
  nTerm=0.75-1.27*fTerm;
  log_10_Pr-=(0.4+0.67*fTerm);                // log10(Pr) + c
  log_10_Pr=log_10_Pr/(nTerm-0.14*log_10_Pr); // d = 0.14
  log_10_Pr*=log_10_Pr;
  fTerm/=(1.0+log_10_Pr);
  fTerm=pow(10.0,fTerm);
  Pcorr = fTerm*Pr/(1.0+Pr);

  K_forward_check = Khigh*Pcorr;
  EXPECT_TRUE(NearScalar(K_forward,K_forward_check)) << "K_forward - three term troe";
  EXPECT_TRUE(NearScalar(K_reverse, 0.0)) << "K_reverse - three term troe";

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


static double GetConcentrationSum(zerork::mechanism *mech,
                              const double temperature,
                              const double pressure,
                              const std::vector<double> &mole_fractions)
{
  const double gas_constant = mech->getGasConstant();
  const int num_species = mech->getNumSpecies();
  const int num_reactions = mech->getNumReactions();

  double Csum = 0;
  for(int j=0; j<num_species; ++j) {
    Csum += mole_fractions[j]*pressure/(gas_constant*temperature);
  }

  return Csum;
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

