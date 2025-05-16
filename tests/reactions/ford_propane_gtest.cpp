#include <math.h>
#include <vector>
#include <cstdlib>

#include <zerork/mechanism.h>

#include <gtest/gtest.h>

// ---------------------------------------------------------------------------
// test constants
// ---------------------------------------------------------------------------
static const double OK_DOUBLE = 1.0e-13; // acceptable relative tolerance

static const char MECH_FILENAME[]  = "mechanisms/real_order/ford_propane.mech";
static const char THERM_FILENAME[] = "mechanisms/ideal/const_specific_heat.therm";
static const char PARSER_LOGNAME[] = "parser.log";

static const int NUM_SPECIES = 6;
static const int NUM_REACTIONS = 2;
static const int NUM_STEPS = 3;

static const int NUM_TEMPS = 5;
static const double TEMP[] = {600.0, 800.0, 1000.0, 1200.0, 2000.0}; 
static const double PRES = 10.0*1.01325e5;  // [Pa]
static const double MOLE_FRAC[] = {0.233, 
                                   0.105, 
                                   0.206, 
                                   0.228, 
                                   0.121, 
                                   0.107};

// ! Two-step propane mechanism inferred from Eqs. (2-5) and Table 4
// !
// ! A. Ghani, T. Poinsot, L. Gicquel, and G. Staffelbach, "LES of longitudinal
// ! and transverse self-excited combustion instabilities in a bluff-body
// ! stabilized turbulent premixed flame," Combustion and Flame, 162(11),
// ! pp. 4075-4083, 2015. DOI:10.1016/j.combustflame.2015.08.024
// ! * downloaded from https://hal.science/hal-01235018/document
//
// ! Arrhenius parameters are given in Table 4 in [cgs] units for the
// ! pre-exponential factor and [cal/mol] for the activation energy, which
// ! is consistent with the default Chemkin units. The stoichiometric
// ! coefficients are given in equations (2-3) and the reactant species
// ! reaction orders are inferred from equations (4-5).
// c3h8 + 3.5 o2  => 3 co + 4 h2o        2.0e12  0.0  3.3e3
// ford /c3h8  0.9028/
// ford /o2    0.6855/
// co + 0.5 o2 <=> co2                   4.5e10  0.0  1.2e3

static double FUEL_ORDER = 0.9028;
static double OXID_ORDER = 0.6855;
static double STEP_ORDER[] = {FUEL_ORDER+OXID_ORDER, 1.5, 1.0};
static double REAC_MOLES[] = {4.5, 1.5, 1.0};
static double PROD_MOLES[] = {7.0, 1.0, 1.5};

static const double LOCAL_GAS_CONST =  8.314462618e3;
static const double EACT_TEMP = 4184.0/LOCAL_GAS_CONST; // convert cal/mol to K
 

static double CalcEquilibriumConstReactionCO_O2(const double temperature);
static void CalcRateCoefficients(const double temperature,
                              std::vector<double> &rate_coef);

static void CalcRateOfProgress(zerork::mechanism *mech,
                              const double temperature,
                              std::vector<double> &concentration,
                              std::vector<double> &rate_of_progress);

static void CalcSpeciesRates(zerork::mechanism *mech,
                              const double temperature,
                              std::vector<double> &concentration,
                              std::vector<double> &spc_create,
                              std::vector<double> &spc_destroy,
                              std::vector<double> &spc_net);

static double GetMechanismRateConstant(zerork::mechanism *mech,
                              const double temperature,
                              const double pressure,
                              const std::vector<double> &mole_fractions,
                              const int reaction_id,
                              const int reaction_dir);

//static void GetUniformComposition(zerork::mechanism *mech,
//                                 std::vector<double> *composition);

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
      load_mech = std::string("../../data/") + load_mech;
      load_therm = std::string("../../data/") + load_therm;
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

TEST_F (MechanismTestFixture, OrderCheck)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";

  zerork::NonIntegerReactionNetwork *net = 
    mechanism_->getNonIntegerReactionNetwork();

  ASSERT_TRUE(net != NULL) <<
    "mechanism_->getNonIntegerReactionNetwork()";

  ASSERT_TRUE(net->GetNumNonIntegerSteps() == NUM_STEPS) <<
    "net->GetNumNonIntegerSteps()";

  for(int j=0; j<NUM_STEPS; ++j) {
     ASSERT_TRUE(NearScalar(net->GetOrderOfStep(j), STEP_ORDER[j])) <<
        "net->GetOrderOfStep(" << j << ")";
  } 

}

TEST_F (MechanismTestFixture, ReactantMoleCheck)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) <<
    "mechanism_ = new zerork::mechanism";

  zerork::NonIntegerReactionNetwork *net =
    mechanism_->getNonIntegerReactionNetwork();

  ASSERT_TRUE(net != NULL) <<
    "mechanism_->getNonIntegerReactionNetwork()";

  ASSERT_TRUE(net->GetNumNonIntegerSteps() == NUM_STEPS) <<
    "net->GetNumNonIntegerSteps()";

  for(int j=0; j<NUM_STEPS; ++j) {
     ASSERT_TRUE(NearScalar(net->GetNumReactantMolesOfStep(j), 
        REAC_MOLES[j])) <<
        "net->GetNumReactantMolesOfStep(" << j << ")";
  } 

}

TEST_F (MechanismTestFixture, ProductMoleCheck)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) <<
    "mechanism_ = new zerork::mechanism";

  zerork::NonIntegerReactionNetwork *net =
    mechanism_->getNonIntegerReactionNetwork();

  ASSERT_TRUE(net != NULL) <<
    "mechanism_->getNonIntegerReactionNetwork()";

  ASSERT_TRUE(net->GetNumNonIntegerSteps() == NUM_STEPS) <<
    "net->GetNumNonIntegerSteps()";

  for(int j=0; j<NUM_STEPS; ++j) {
     ASSERT_TRUE(NearScalar(net->GetNumProductMolesOfStep(j),
        PROD_MOLES[j])) <<
        "net->GetNumProductMolesOfStep(" << j << ")";
  } 

}

TEST_F (MechanismTestFixture, RateCoef)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) <<
    "mechanism_ = new zerork::mechanism";

  std::vector<double> calc_rate_coef, mech_fwd_coef, mech_rev_coef;
  std::vector<double> concentration;

  mech_fwd_coef.assign(NUM_REACTIONS, 0.0);
  mech_rev_coef.assign(NUM_REACTIONS, 0.0);
  concentration.assign(NUM_SPECIES, 0.0);

  for(int j=0; j<NUM_TEMPS; ++j) {

     double temperature = TEMP[j];
     double total_conc = PRES/(LOCAL_GAS_CONST*temperature);
     
     for(int k=0; k<NUM_SPECIES; ++k) {
       concentration[k] = MOLE_FRAC[k] * total_conc;
     }

     //printf("T = %.0f [K]\n",temperature);

     CalcRateCoefficients(temperature, calc_rate_coef);
     
     mechanism_->getKrxnFromTC(temperature, &concentration[0],
       &mech_fwd_coef[0],
       &mech_rev_coef[0]);

     ASSERT_TRUE(NearScalar(calc_rate_coef[0], mech_fwd_coef[0])) <<
       "reaction 0 (fwd)" << " at T = " << temperature << " [K]";

     ASSERT_TRUE(NearScalar(calc_rate_coef[1], mech_fwd_coef[1])) <<
       "reaction 1 (fwd)" << " at T = " << temperature << " [K]";

     ASSERT_TRUE(NearScalar(calc_rate_coef[2], mech_rev_coef[1])) <<
       "reaction 1 (rev)" << " at T = " << temperature << " [K]";
  }
}

TEST_F (MechanismTestFixture, RateOfProgress)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) <<
    "mechanism_ = new zerork::mechanism";

  std::vector<double> calc_step_rop, mech_net_rate, mech_create_rate;
  std::vector<double> mech_destroy_rate, mech_step_rop;

  std::vector<double> concentration;


  mech_net_rate.assign(NUM_SPECIES, 0.0);
  mech_create_rate.assign(NUM_SPECIES, 0.0);
  mech_destroy_rate.assign(NUM_SPECIES, 0.0);
  mech_step_rop.assign(NUM_STEPS, 0.0);
 
  concentration.assign(NUM_SPECIES, 0.0);

  for(int j=0; j<NUM_TEMPS; ++j) {

     double temperature = TEMP[j];
     double total_conc = PRES/(LOCAL_GAS_CONST*temperature);

     for(int k=0; k<NUM_SPECIES; ++k) {
       concentration[k] = MOLE_FRAC[k] * total_conc;
     }

     //printf("T = %.0f [K]\n",temperature);

     CalcRateOfProgress(mechanism_,temperature, concentration, calc_step_rop);

     mechanism_->getReactionRates(temperature, &concentration[0],
       &mech_net_rate[0],
       &mech_create_rate[0],
       &mech_destroy_rate[0],
       &mech_step_rop[0]);

     for(int k=0; k<NUM_STEPS; ++k) {

        ASSERT_TRUE(NearScalar(calc_step_rop[k], mech_step_rop[k])) <<
          "step " << k << " at T = " << temperature << " [K]";
     }
  }
}

TEST_F (MechanismTestFixture, SpeciesRates)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) <<
    "mechanism_ = new zerork::mechanism";

  std::vector<double> calc_net_rate, calc_create_rate, calc_destroy_rate; 
  std::vector<double> mech_net_rate, mech_create_rate;
  std::vector<double> mech_destroy_rate, mech_step_rop;

  std::vector<double> concentration;

  calc_net_rate.assign(NUM_SPECIES, 0.0);
  calc_create_rate.assign(NUM_SPECIES, 0.0);
  calc_destroy_rate.assign(NUM_SPECIES, 0.0);
  mech_net_rate.assign(NUM_SPECIES, 0.0);
  mech_create_rate.assign(NUM_SPECIES, 0.0);
  mech_destroy_rate.assign(NUM_SPECIES, 0.0);

  mech_step_rop.assign(NUM_STEPS, 0.0);

  concentration.assign(NUM_SPECIES, 0.0);

  for(int j=0; j<NUM_TEMPS; ++j) {

     double temperature = TEMP[j];
     double total_conc = PRES/(LOCAL_GAS_CONST*temperature);

     for(int k=0; k<NUM_SPECIES; ++k) {
       concentration[k] = MOLE_FRAC[k] * total_conc;
     }

     //printf("T = %.0f [K]\n",temperature);

     CalcSpeciesRates(mechanism_, temperature, concentration,
                      calc_create_rate,
                      calc_destroy_rate,
                      calc_net_rate);

     mechanism_->getReactionRates(temperature, &concentration[0],
       &mech_net_rate[0],
       &mech_create_rate[0],
       &mech_destroy_rate[0],
       &mech_step_rop[0]);

     for(int k=0; k<NUM_SPECIES; ++k) {

        ASSERT_TRUE(NearScalar(calc_create_rate[k], mech_create_rate[k])) <<
          "step " << k << " at T = " << temperature << " [K]";
        ASSERT_TRUE(NearScalar(calc_destroy_rate[k], mech_destroy_rate[k])) <<
          "step " << k << " at T = " << temperature << " [K]";
        ASSERT_TRUE(NearScalar(calc_net_rate[k], mech_net_rate[k])) <<
          "step " << k << " at T = " << temperature << " [K]";

     }
  }
}




int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

static double CalcEquilibriumConstReactionCO_O2(const double temperature)
{
  // reaction 1: co + 0.5 o2 <=> co2
  // S/Ru = a_0*log(T) + a_6
  double co_entropy   = ( 3.50000000e+00)*log(temperature) + 3.83145325e+00;
  double o2_entropy   = ( 3.50000000e+00)*log(temperature) + 4.73253404e+00;
  double co2_entropy  = ( 4.46541338e+00)*log(temperature) + 2.70309557e-01;
  // H/RuT = a_0 + a_5/T
  double co_enthalpy  = 3.50000000e+00 - 1.43372329e+04/temperature;
  double o2_enthalpy  = 3.50000000e+00 - 1.04352500e+03/temperature;
  double co2_enthalpy = 4.46541338e+00 - 4.86597535e+04/temperature;

  double delta_entropy  = -co_entropy  - 0.5*o2_entropy  + co2_entropy;
  double delta_enthalpy = -co_enthalpy - 0.5*o2_enthalpy + co2_enthalpy;

  double Keq = pow(1.01325e5/(LOCAL_GAS_CONST*temperature), -0.5)*
    exp(delta_entropy-delta_enthalpy);

  return Keq;
}

static void CalcRateCoefficients(const double temperature,
                              std::vector<double> &rate_coef)
{
   rate_coef.clear();
   rate_coef.assign(NUM_STEPS, 0.0);

   // step 0: c3h8 + 3.5 o2  => 3 co + 4 h2o        2.0e12  0.0  3.3e3
   rate_coef[0] = 2.0e12  *  exp(-EACT_TEMP *   3.3e3  /temperature); 
   // units [(mol/cm**3)**N /s], with N = -0.5883 = 1 - (0.9028+0.6855)

   // step 1: co + 0.5 o2 <=> co2                   4.5e10  0.0  1.2e3
   rate_coef[1] = 4.5e10  *  exp(-EACT_TEMP *   1.2e3  /temperature);
   // units [(mol/cm**3)**N /s], with N = -0.5 = 1 - (1+0.5)

   // zero-rk units [kmol/m**3]**N
   // 1 [mol/cm**3] = 1 * (0.001 kmol/mol) * (100 cm/m)**3
   // 1 [mol/cm**3] = 1000 [kmol/m**3]
   
   // convert step 0 and 1 to zero-rk units
   rate_coef[0] *= pow(1000.0, 1.0-STEP_ORDER[0]);
   rate_coef[1] *= pow(1000.0, 1.0-STEP_ORDER[1]);

   // calculate step 2 using equilibrium K_rev = K_fwd/K_eqm
   double K_eqm = CalcEquilibriumConstReactionCO_O2(temperature);
   rate_coef[2] = rate_coef[1]/K_eqm;

}

static void CalcRateOfProgress(zerork::mechanism *mech,
                              const double temperature,
                              std::vector<double> &concentration,
                              std::vector<double> &rate_of_progress)
{
  std::vector<double> rate_coef;

  rate_of_progress.clear();
  rate_of_progress.assign(NUM_STEPS, 0.0);
  
  CalcRateCoefficients(temperature, rate_coef);

  
  double C_c3h8 = concentration[mech->getIdxFromName("c3h8")];
  double C_o2   = concentration[mech->getIdxFromName("o2")];
  double C_co   = concentration[mech->getIdxFromName("co")];
  double C_co2  = concentration[mech->getIdxFromName("co2")];
  double C_h2o  = concentration[mech->getIdxFromName("h2o")];

  // c3h8 + 3.5 o2  => 3 co + 4 h2o        2.0e12  0.0  3.3e3
  // ford /c3h8  0.9028/
  // ford /o2    0.6855/
  // co + 0.5 o2 <=> co2                   4.5e10  0.0  1.2e3

  rate_of_progress[0] = rate_coef[0] * pow(C_c3h8, FUEL_ORDER)
    * pow(C_o2, OXID_ORDER);

  rate_of_progress[1] = rate_coef[1] * C_co * pow(C_o2, 0.5);

  rate_of_progress[2] = rate_coef[2] * C_co2; 

}

static void CalcSpeciesRates(zerork::mechanism *mech,
                              const double temperature,
                              std::vector<double> &concentration,
                              std::vector<double> &spc_create,
                              std::vector<double> &spc_destroy,
                              std::vector<double> &spc_net)
{
  std::vector<double> rate_of_progress;
  
  const int id_n2   = mech->getIdxFromName("n2");
  const int id_o2   = mech->getIdxFromName("o2");
  const int id_c3h8 = mech->getIdxFromName("c3h8");
  const int id_co   = mech->getIdxFromName("co");
  const int id_co2  = mech->getIdxFromName("co2");
  const int id_h2o  = mech->getIdxFromName("h2o");

  spc_create.clear();
  spc_destroy.clear();
  spc_net.clear();

  spc_create.assign(NUM_SPECIES, 0.0);
  spc_destroy.assign(NUM_SPECIES, 0.0);
  spc_net.assign(NUM_SPECIES, 0.0);

 
  CalcRateOfProgress(mech, temperature, concentration, rate_of_progress);

  // c3h8 + 3.5 o2  => 3 co + 4 h2o        2.0e12  0.0  3.3e3
  // ford /c3h8  0.9028/
  // ford /o2    0.6855/
  // co + 0.5 o2 <=> co2                   4.5e10  0.0  1.2e3

  // creation rates
  spc_create[id_co]    += 3.0 * rate_of_progress[0]; // reaction 0 (fwd)
  spc_create[id_h2o]   += 4.0 * rate_of_progress[0];

  spc_create[id_co2]   += 1.0 * rate_of_progress[1]; // reaction 1 (fwd)

  spc_create[id_co]    += 1.0 * rate_of_progress[2]; // reaction 1 (rev)
  spc_create[id_o2]    += 0.5 * rate_of_progress[2];

  // destruction rates
  spc_destroy[id_c3h8] += 1.0 * rate_of_progress[0]; // reaction 0 (fwd)
  spc_destroy[id_o2]   += 3.5 * rate_of_progress[0]; 

  spc_destroy[id_co]   += 1.0 * rate_of_progress[1]; // reaction 1 (fwd)
  spc_destroy[id_o2]   += 0.5 * rate_of_progress[1];

  spc_destroy[id_co2]  += 1.0 * rate_of_progress[2]; // reaction 1 (rev)

  // net rates
  for(int j=0; j<NUM_SPECIES; ++j) {
    spc_net[j] = spc_create[j] - spc_destroy[j];
  }

}



static double GetMechanismRateConstant(zerork::mechanism *mech,
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

// static double GetEquilibriumConst(const double temperature,
//                                   const double pressure)
// {
//   // A + A <=> A2 (diatomic)
//   // S/Ru
//   double   atomic_entropy =   ATOMIC_A0*log(temperature) +   ATOMIC_A6;
//   double diatomic_entropy = DIATOMIC_A0*log(temperature) + DIATOMIC_A6;
//   // H/RuT
//   double   atomic_enthalpy =   ATOMIC_A0 +   ATOMIC_A5/temperature;
//   double diatomic_enthalpy = DIATOMIC_A0 + DIATOMIC_A5/temperature;

//   double delta_entropy  = -2.0*atomic_entropy  + diatomic_entropy;
//   double delta_enthalpy = -2.0*atomic_enthalpy + diatomic_enthalpy;

//   double Keq = pow(1.01325e5/(LOCAL_GAS_CONST*temperature), -1.0)*
//     exp(delta_entropy-delta_enthalpy);

//   return Keq;
// }





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

// static bool NearVector(size_t n,
//                        const double a[],
//                        const double b[])
// {
//   bool near_scalar;
//   for(size_t j=0; j<n; ++j) {
//     near_scalar = NearScalar(a[j],b[j]);
//     if(near_scalar == false) {
//       return false;
//     }
//   }
//   return true;
// }
