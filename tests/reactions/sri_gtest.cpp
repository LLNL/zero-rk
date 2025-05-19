#include <math.h>
#include <vector>
#include <cstdlib>

#include <zerork/mechanism.h>

#include <gtest/gtest.h>

// ---------------------------------------------------------------------------
// test constants
// ---------------------------------------------------------------------------
static const double OK_DOUBLE = 5.0e-15; // acceptable relative tolerance

static const char MECH_FILENAME[]  = "mechanisms/ideal/sri_reaction.mech";
static const char THERM_FILENAME[] = "mechanisms/ideal/const_specific_heat.therm";
static const char PARSER_LOGNAME[] = "parser.log";

static const int NUM_SPECIES = 3;
static const int NUM_REACTIONS = 2;
static const int NUM_STEPS = 4;


static const double LOCAL_GAS_CONST = 8314.462618;

static double GetRateConstant(zerork::mechanism *mech,
                              const double temperature,
                              const double pressure,
                              const std::vector<double> &mole_fractions,
                              const int reaction_id,
                              const int reaction_dir);

static void GetUniformComposition(zerork::mechanism *mech,
                                  std::vector<double> *composition);

static double KForwardReaction01(const double temperature,
                                 const double pressure);

static double KForwardReaction02(const double temperature,
                                 const double pressure);


static bool NearScalar(const double a,
                      const double b);
//static bool NearVector(size_t n,
//                       const double a[],
//                      const double b[]);

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
  const int num_steps     = mechanism_->getNumSteps();

  ASSERT_TRUE(NUM_SPECIES == num_species) <<
    "mechanism_->getNumSpecies()";

  ASSERT_TRUE(NUM_REACTIONS == num_reactions) <<
    "mechanism_->getNumReactions()";

  ASSERT_TRUE(NUM_STEPS == num_steps) <<
    "mechanism_->getNumSteps()";
}

TEST_F (MechanismTestFixture, Reaction01_Three_Param_SRI)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";

  const double p_mult = 2.0;
  const double delta_t = 100.0;
  const int num_pressures = 24;
  const int num_temperatures = 36;

  double initial_pressure = 1.0;
  double temperature = 300.0; 
  double pressure;
  std::vector<double> mole_fractions;
  double k_forward_zero_rk;

  GetUniformComposition(mechanism_,
                        &mole_fractions);

  for(int j=0; j<num_temperatures; ++j) {
    pressure = initial_pressure;
    for(int k=0; k<num_pressures; ++k) {

      k_forward_zero_rk = GetRateConstant(mechanism_,
                                          temperature,
                                          pressure,
                                          mole_fractions,
                                          0,  // reaction index
                                          1); // forward direction
      EXPECT_TRUE(NearScalar(k_forward_zero_rk,
                             KForwardReaction01(temperature,pressure))) <<
        "Reaction 01 forward rate constant at T = " << temperature <<
        " [K] and p = " << pressure << " [Pa]" << endl;

      pressure *= p_mult;
    }
    temperature += delta_t;
  }
  
}

TEST_F (MechanismTestFixture, Reaction02_Five_Param_SRI)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";

  const double p_mult = 2.0;
  const double delta_t = 100.0;
  const int num_pressures = 24;
  const int num_temperatures = 36;

  double initial_pressure = 1.0;
  double temperature = 300.0; 
  double pressure;
  std::vector<double> mole_fractions;
  double k_forward_zero_rk;

  GetUniformComposition(mechanism_,
                        &mole_fractions);

  for(int j=0; j<num_temperatures; ++j) {
    pressure = initial_pressure;
    for(int k=0; k<num_pressures; ++k) {

      k_forward_zero_rk = GetRateConstant(mechanism_,
                                          temperature,
                                          pressure,
                                          mole_fractions,
                                          1,  // reaction index
                                          1); // forward direction
      EXPECT_TRUE(NearScalar(k_forward_zero_rk,
                             KForwardReaction02(temperature,pressure))) <<
        "Reaction 02 forward rate constant at T = " << temperature <<
        " [K] and p = " << pressure << " [Pa]" << endl;

      pressure *= p_mult;
    }
    temperature += delta_t;
  }
  
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
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

static double KForwardReaction01(const double temperature,
                                 const double pressure) 
{
  // the factor 1 kmol/m^3 * 1.0e-6 m^3/cm^3 * 1000.0 mol/kmol = 0.001 mol/cm^3
  double M = 0.001*pressure/(LOCAL_GAS_CONST*temperature); // [mol/cm^3]
  double k_low = 8.0e26*pow(temperature, -3.0);  // [cm^6/mol^2/s]
  double k_inf = 6.0e16*pow(temperature, -1.0);  // [cm^3/mol/s]
  double Pr = k_low*M/k_inf;  // dimensionless
  double a = 0.45;
  double b = 797.0; // [K]
  double c = 979.0; // [K]
  double d = 1.0;
  double e = 0.0;
  double log_10_Pr = log10(Pr);
  double x_power = 1.0/(1.0+log_10_Pr*log_10_Pr);

  double F = d*pow(a*exp(-b/temperature) + exp(-temperature/c),x_power)*
    pow(temperature, e);
  // the factor 1 cm^3/mol * 1.0e-6 m^3/cm^3 * 1000.0 mol/kmol = 0.001 m^3/kmol
  return k_inf*0.001*(Pr/(1.0+Pr))*F;
}

static double KForwardReaction02(const double temperature,
                                 const double pressure) 
{
  // the factor 1 kmol/m^3 * 1.0e-6 m^3/cm^3 * 1000.0 mol/kmol = 0.001 mol/cm^3
  double M = 0.001*pressure/(LOCAL_GAS_CONST*temperature);
  double k_low = 8.0e26*pow(temperature, -3.0);  // [cm^6/mol^2/s]
  double k_inf = 6.0e16*pow(temperature, -1.0);  // [cm^3/mol/s]
  double Pr = k_low*M/k_inf;  // dimensionless
  double a = 0.45;
  double b = 797.0; // [K]
  double c = 979.0; // [K]
  double d = 0.9;
  double e = 2.0;
  double log_10_Pr = log10(Pr);
  double x_power = 1.0/(1.0+log_10_Pr*log_10_Pr);

  double F = d*pow(a*exp(-b/temperature) + exp(-temperature/c),x_power)*
    pow(temperature, e);
  // the factor 1 cm^3/mol = 1.0e-6 m^3/cm^3 * 1000.0 mol/kmol = 0.001 m^3/kmol
  return k_inf*0.001*(Pr/(1.0+Pr))*F;
}
