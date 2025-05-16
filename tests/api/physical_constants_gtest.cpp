#include <math.h>

#include <gtest/gtest.h>

#include <zerork/constants_api.h>

// ---------------------------------------------------------------------------
// test constants
// ---------------------------------------------------------------------------
static const double OK_DOUBLE = 1.0e-3; // acceptable relative tolerance

static bool NearScalar(const double a,
                       const double b,
                       const double rel_tol,
                       const double abs_tol);

// ---------------------------------------------------------------------------
// test fixture for the ZeroRK PhysicalConstants class
class PhysicalConstantsTestFixture: public ::testing::Test 
{ 
 public: 
  PhysicalConstantsTestFixture( ) { 
    // initialization code here
    physical_constants_  = NULL;
    physical_constants_ = new zerork::PhysicalConstants();
  } 

  void SetUp( ) { 
    // code here will execute just before the test ensues 
  }

  void TearDown( ) { 
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be
  }

  ~PhysicalConstantsTestFixture( )  { 
    // cleanup any pending stuff, but no exceptions allowed
    if(physical_constants_ != NULL) {
      delete physical_constants_;
    }
  }

  // put in any custom data members that you need 
  zerork::PhysicalConstants *physical_constants_;
};

TEST_F (PhysicalConstantsTestFixture, Allocation)
{
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(physical_constants_ != NULL) << 
    "physical_constants_ = new PhysicalConstants()";

}

TEST_F (PhysicalConstantsTestFixture, GasConstant)
{
  double relative_tolerance = OK_DOUBLE;
  double absolute_tolerance = 1.0e-20;
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(physical_constants_ != NULL) << 
    "physical_constants_ = new PhysicalConstants()";

  EXPECT_TRUE(NearScalar(physical_constants_->GetGasConstant(),
                         8.314462618e3,
                         relative_tolerance,
                         absolute_tolerance)) << std::setprecision(3) <<
    "with relative_tolerance = "  << relative_tolerance <<
    ", and absolute_tolerance = " << absolute_tolerance;

}

TEST_F (PhysicalConstantsTestFixture, AvogadroNum)
{
  double relative_tolerance = OK_DOUBLE;
  double absolute_tolerance = 1.0e-20;
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(physical_constants_ != NULL) << 
    "physical_constants_ = new PhysicalConstants()";

  EXPECT_TRUE(NearScalar(physical_constants_->GetAvogadroNum(),
                         6.022140767e26,
                         relative_tolerance,
                         absolute_tolerance)) << std::setprecision(3) <<
    "with relative_tolerance = "  << relative_tolerance <<
    ", and absolute_tolerance = " << absolute_tolerance;

}
TEST_F (PhysicalConstantsTestFixture, BoltzmannConstant)
{
  double relative_tolerance = OK_DOUBLE;
  double absolute_tolerance = 1.0e-20;
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(physical_constants_ != NULL) << 
    "physical_constants_ = new PhysicalConstants()";

  EXPECT_TRUE(NearScalar(physical_constants_->GetBoltzmannConstant(),
                         1.380649e-23,
                         relative_tolerance,
                         absolute_tolerance)) << std::setprecision(3) <<
    "with relative_tolerance = "  << relative_tolerance <<
    ", and absolute_tolerance = " << absolute_tolerance;

}
TEST_F (PhysicalConstantsTestFixture, JoulesPerCalorie)
{
  double relative_tolerance = OK_DOUBLE;
  double absolute_tolerance = 1.0e-20;
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(physical_constants_ != NULL) << 
    "physical_constants_ = new PhysicalConstants()";

  EXPECT_TRUE(NearScalar(physical_constants_->GetJoulesPerCalorie(),
                         4184.0,
                         relative_tolerance,
                         absolute_tolerance)) << std::setprecision(3) <<
    "with relative_tolerance = "  << relative_tolerance <<
    ", and absolute_tolerance = " << absolute_tolerance;

}
TEST_F (PhysicalConstantsTestFixture, HydrogenMass)
{
  double relative_tolerance = OK_DOUBLE;
  double absolute_tolerance = 1.0e-20;
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(physical_constants_ != NULL) << 
    "physical_constants_ = new PhysicalConstants()";

  EXPECT_TRUE(NearScalar(physical_constants_->GetAtomicMass(1),
                         1.00794,
                         relative_tolerance,
                         absolute_tolerance)) << std::setprecision(3) <<
    "with relative_tolerance = "  << relative_tolerance <<
    ", and absolute_tolerance = " << absolute_tolerance;

}
TEST_F (PhysicalConstantsTestFixture, CarbonMass)
{
  double relative_tolerance = OK_DOUBLE;
  double absolute_tolerance = 1.0e-20;
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(physical_constants_ != NULL) << 
    "physical_constants_ = new PhysicalConstants()";

  EXPECT_TRUE(NearScalar(physical_constants_->GetAtomicMass(6),
                         12.0107,
                         relative_tolerance,
                         absolute_tolerance)) << std::setprecision(3) <<
    "with relative_tolerance = "  << relative_tolerance <<
    ", and absolute_tolerance = " << absolute_tolerance;

}
TEST_F (PhysicalConstantsTestFixture, NitrogenMass)
{
  double relative_tolerance = OK_DOUBLE;
  double absolute_tolerance = 1.0e-20;
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(physical_constants_ != NULL) << 
    "physical_constants_ = new PhysicalConstants()";

  EXPECT_TRUE(NearScalar(physical_constants_->GetAtomicMass(7),
                         14.0067,
                         relative_tolerance,
                         absolute_tolerance)) << std::setprecision(3) <<
    "with relative_tolerance = "  << relative_tolerance <<
    ", and absolute_tolerance = " << absolute_tolerance;

}
TEST_F (PhysicalConstantsTestFixture, OxygenMass)
{
  double relative_tolerance = OK_DOUBLE;
  double absolute_tolerance = 1.0e-20;
  // Note assert failures halts test fixture, but not all tests
  ASSERT_TRUE(physical_constants_ != NULL) << 
    "physical_constants_ = new PhysicalConstants()";

  EXPECT_TRUE(NearScalar(physical_constants_->GetAtomicMass(8),
                         15.9994,
                         relative_tolerance,
                         absolute_tolerance)) << std::setprecision(3) <<
    "with relative_tolerance = "  << relative_tolerance <<
    ", and absolute_tolerance = " << absolute_tolerance;

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
