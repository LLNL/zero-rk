#include <math.h>

#include <gtest/gtest.h>

#include <zerork/mechanism.h>

// ---------------------------------------------------------------------------
// test constants
// ---------------------------------------------------------------------------
static const double OK_DOUBLE = 1.0e-3; // acceptable relative tolerance

static const char MECH_FILENAME[]  = "mechanisms/ideal/plog_test.mech";
static const char THERM_FILENAME[] = "mechanisms/ideal/const_specific_heat.therm";
static const char PARSER_LOGNAME[] = "parser.log";

static bool NearScalar(const double a,
                       const double b,
                       const double rel_tol,
                       const double abs_tol);

// ---------------------------------------------------------------------------
// test fixture for the ZeroRK mechanism class
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
  zerork::mechanism *mechanism_;
};

TEST_F (MechanismTestFixture, Allocation)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) <<
    "mechanism_ = new mechanism()";
}


TEST_F (MechanismTestFixture, ElementInfo)
{
  double relative_tolerance = OK_DOUBLE;
  double absolute_tolerance = 1.0e-20;
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) <<
    "mechanism_ = new mechanism()";

  std::map<std::string, double> elementInfo = mechanism_->getElementInfo();
  EXPECT_TRUE(elementInfo.size() == 1);
  EXPECT_TRUE(elementInfo.begin()->first == std::string("H"));
  EXPECT_TRUE(NearScalar(elementInfo.begin()->second,
                         1.00794,
                         relative_tolerance,
                         absolute_tolerance)) << std::setprecision(3) <<
    "with relative_tolerance = "  << relative_tolerance <<
    ", and absolute_tolerance = " << absolute_tolerance;

}

TEST_F (MechanismTestFixture, SpeciesElementInfo)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(mechanism_ != NULL) <<
    "mechanism_ = new mechanism()";

  std::map<std::string, std::map<std::string, int> > speciesElementInfo = mechanism_->getSpeciesElementInfo();
  EXPECT_TRUE(speciesElementInfo.size() == 2);
  EXPECT_TRUE(speciesElementInfo.begin()->first == std::string("h"));
  EXPECT_TRUE(speciesElementInfo.begin()->second.begin()->first== std::string("H"));
  EXPECT_TRUE(speciesElementInfo.begin()->second.begin()->second== 1);
  EXPECT_TRUE((++speciesElementInfo.begin())->first == std::string("h2"));
  EXPECT_TRUE((++speciesElementInfo.begin())->second.begin()->first== std::string("H"));
  EXPECT_TRUE((++speciesElementInfo.begin())->second.begin()->second== 2);

}


// --------------------------------------------------------------------------

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

static bool NearScalar(const double a,
                       const double b,
                       const double rel_tol,
                       const double abs_tol)
{
  double weight = 0.5*(fabs(a)+fabs(b));

  if(weight > fabs(abs_tol)) {
    // check the normalized difference
    if(fabs(a-b)/weight > fabs(rel_tol)) {

      printf("# NearScalar false: %24.18e != %24.18e\n",a,b);
      return false;
    }
  }
  return true;
}
