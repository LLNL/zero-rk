#include <math.h>
#include <vector>
#include <cstdlib>

#include <zerork/mechanism.h>

#include <gtest/gtest.h>

// ---------------------------------------------------------------------------
// test constants
// ---------------------------------------------------------------------------
static const double OK_DOUBLE = 1.0e-15; // acceptable relative tolerance

static const char MECH_FILENAME[]  = "mechanisms/ideal/big_molecule.mech";
static const char THERM_FILENAME[] = "mechanisms/ideal/big_molecule.therm";
static const char PARSER_LOGNAME[] = "parser.log";

static const int NUM_SPECIES = 13;
static const int NUM_REACTIONS = 4;
static const int NUM_STEPS = 4;

// Masses [kg/kmol] obtained from the publicshed IUPAC 2007 values
static const double H_WEIGHT =  1.00794;
static const double C_WEIGHT = 12.0107;
static const double O_WEIGHT = 15.9994;

static const double LOCAL_GAS_CONST = 8314.4598;

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
TEST_F (MechanismTestFixture, LargeMoleculeMolecularWeights)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) << 
    "mechanism_ = new zerork::mechanism";

  std::vector<double>molecular_weights;
  molecular_weights.assign(NUM_SPECIES, 0.0);

  mechanism_->getMolWtSpc(&molecular_weights[0]);
  
  EXPECT_TRUE(NearScalar(molecular_weights[0],
                         H_WEIGHT)) << "Molecular weight of H";

  EXPECT_TRUE(NearScalar(molecular_weights[1],
                         O_WEIGHT+H_WEIGHT)) << "Molecular weight of OH";

  EXPECT_TRUE(NearScalar(molecular_weights[2],
                         2*C_WEIGHT+2*H_WEIGHT)) << "Molecular weight of C_2H_2";

  EXPECT_TRUE(NearScalar(molecular_weights[3],
                         2.5e6*C_WEIGHT+4.99999e5*H_WEIGHT)) << "Molecular weight of BIN18AJ (C_2500000H_499999)";

  EXPECT_TRUE(NearScalar(molecular_weights[4],
                         2.5e6*C_WEIGHT+5.0e5*H_WEIGHT)) << "Molecular weight of BIN18A (C_2500000H_500000)";

  EXPECT_TRUE(NearScalar(molecular_weights[5],
                         2.5e6*C_WEIGHT+1.24999e5*H_WEIGHT)) << "Molecular weight of BIN18BJ (C_2500000H_124999)";

  EXPECT_TRUE(NearScalar(molecular_weights[6],
                         2.5e6*C_WEIGHT+1.25e5*H_WEIGHT)) << "Molecular weight of BIN18B (C_2500000H_125000)";

  EXPECT_TRUE(NearScalar(molecular_weights[7],
                         5.0e6*C_WEIGHT+9.99999e5*H_WEIGHT)) << "Molecular weight of BIN19AJ (C_5000000H_999999)";

  EXPECT_TRUE(NearScalar(molecular_weights[8],
                         5.0e6*C_WEIGHT+1.0e6*H_WEIGHT)) << "Molecular weight of BIN19A (C_5000000H_1000000)";
 
  EXPECT_TRUE(NearScalar(molecular_weights[9],
                         5.0e6*C_WEIGHT+2.49999e5*H_WEIGHT)) << "Molecular weight of BIN19BJ (C_5000000H_249999)";

  EXPECT_TRUE(NearScalar(molecular_weights[10],
                         5.0e6*C_WEIGHT+2.5e5*H_WEIGHT)) << "Molecular weight of BIN19B (C_5000000H_250000)";
 
  EXPECT_TRUE(NearScalar(molecular_weights[11],
                         1.0e7*C_WEIGHT+2.0e6*H_WEIGHT)) << "Molecular weight of BIN20A (C_10000000H_249999)";

  EXPECT_TRUE(NearScalar(molecular_weights[12],
                         1.0e7*C_WEIGHT+5.0e5*H_WEIGHT)) << "Molecular weight of BIN20B (C_10000000H_250000)";
  

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
