#include <math.h>

#include <gtest/gtest.h>

#include <zerork/elemental_composition.h>

// ---------------------------------------------------------------------------
// test constants
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// test fixture for the ZeroRK ElementalComposition class
class ElementalCompositionTestFixture: public ::testing::Test 
{ 
 public: 
  ElementalCompositionTestFixture( )
  :
   test_name1("H2"),
   test_name2("O2"),
   test_name3("C2H30")
  { 
    // initialization code here
  } 

  void SetUp( ) { 
    // code here will execute just before the test ensues 
  }

  void TearDown( ) { 
    // code here will be called just after the test completes
    // ok to through exceptions from here if need be
  }

  ~ElementalCompositionTestFixture( )  { 
  }

  // put in any custom data members that you need 
  std::string test_name1;
  std::string test_name2;
  std::string test_name3;
};

TEST_F (ElementalCompositionTestFixture, InstantiationDefault)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition();
  SUCCEED();
}

TEST_F (ElementalCompositionTestFixture, InstantiationWithName)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition(test_name1);
  SUCCEED();
}

TEST_F (ElementalCompositionTestFixture, InstantiationByCopy)
{
  zerork::ElementalComposition elemental_composition1 = zerork::ElementalComposition();
  zerork::ElementalComposition elemental_composition2 = zerork::ElementalComposition(elemental_composition1);
  EXPECT_NE(&elemental_composition2,  &elemental_composition1) << 
    "&elemental_composition2 != &elemental_composition1";
}

TEST_F (ElementalCompositionTestFixture, CopyOperator)
{
  zerork::ElementalComposition elemental_composition1 = zerork::ElementalComposition();
  zerork::ElementalComposition elemental_composition2 = elemental_composition1;
  EXPECT_NE(&elemental_composition2, &elemental_composition1) << 
    "&elemental_composition2 != &elemental_composition1";
}


TEST_F (ElementalCompositionTestFixture, GetNumAtoms)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition(test_name1);
  elemental_composition["H"] = 2;
  elemental_composition["C"] = 3;
  EXPECT_EQ(elemental_composition.GetNumAtoms("H"), 2) << 
    "elemental_composition.GetNumAtoms(\"H\")";
  EXPECT_EQ(elemental_composition.GetNumAtoms("C"), 3) << 
    "elemental_composition.GetNumAtoms(\"C\")";
}

TEST_F (ElementalCompositionTestFixture, GetNumHeavyAtoms)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition(test_name1);
  elemental_composition["H"] = 2;
  elemental_composition["C"] = 3;
  elemental_composition["O"] = 1;
  EXPECT_EQ(elemental_composition.GetNumHeavyAtoms(), 4) << 
    "elemental_composition.GetNumHeavyAtoms()";
}

TEST_F (ElementalCompositionTestFixture, GetNumTotalAtoms)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition(test_name1);
  elemental_composition["H"] = 2;
  elemental_composition["C"] = 3;
  elemental_composition["O"] = 1;
  EXPECT_EQ(elemental_composition.GetNumTotalAtoms(), 6) << 
    "elemental_composition.GetNumTotalAtoms()";
}

TEST_F (ElementalCompositionTestFixture, GetElementVector)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition(test_name1);
  elemental_composition["H"] = 2;
  elemental_composition["C"] = 3;
  elemental_composition["O"] = 1;
  std::vector<std::string> elements_check;
  elements_check.push_back("C");
  elements_check.push_back("H");
  elements_check.push_back("O");
  EXPECT_EQ(elemental_composition.GetElementVector(), elements_check) << 
    "elemental_composition.GetElementVector()";
}

TEST_F (ElementalCompositionTestFixture, name)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition(test_name1);
  elemental_composition.name() = test_name2;
  EXPECT_EQ(elemental_composition.name(), test_name2) << 
    "elemental_composition.name()";
}

TEST_F (ElementalCompositionTestFixture, name_const)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition(test_name1);
  const std::string name_check = elemental_composition.name();
  EXPECT_EQ(name_check, test_name1) << 
    "elemental_composition.name() const failed";
}

TEST_F (ElementalCompositionTestFixture, ToString)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition();
  elemental_composition["H"] = 2;
  elemental_composition["C"] = 3;
  elemental_composition["O"] = 1;
  std::string string_check("C3H2O");
  EXPECT_EQ(elemental_composition.ToString(), string_check) << 
    "elemental_composition.ToString() produced wrong result";
}

TEST_F (ElementalCompositionTestFixture, ToStringWithSeparator)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition();
  elemental_composition["H"] = 2;
  elemental_composition["C"] = 3;
  elemental_composition["O"] = 1;
  std::string string_check("C=3 H=2 O=1");
  EXPECT_EQ(elemental_composition.ToStringWithSeparator("="), string_check) << 
    "elemental_composition.ToStringWithSeparator(\"=\") produced wrong result.";
}

TEST_F (ElementalCompositionTestFixture, SubscriptOperator)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition();
  elemental_composition["H"] = 2;
  EXPECT_EQ(elemental_composition["H"], 2) << 
    "elemental_composition[\"H\"] != 2";
}

TEST_F (ElementalCompositionTestFixture, SubscriptOperatorConst)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition();
  elemental_composition["H"] = 2;
  const int& check_val = elemental_composition["H"];
  EXPECT_EQ(check_val, 2) << 
    "const int& check_val = elemental_composition[\"H\"]";
}

TEST_F (ElementalCompositionTestFixture, clear)
{
  zerork::ElementalComposition elemental_composition = zerork::ElementalComposition();
  elemental_composition["H"] = 2;
  elemental_composition["C"] = 3;
  elemental_composition["O"] = 1;
  elemental_composition.clear();
  EXPECT_EQ(elemental_composition.GetNumTotalAtoms(), 0) << 
    "elemental_composition.clear()";
}

TEST_F (ElementalCompositionTestFixture, RelationalOperatorsEqual)
{
  zerork::ElementalComposition elemental_composition1 = zerork::ElementalComposition("AA");
  zerork::ElementalComposition elemental_composition2 = elemental_composition1;
  EXPECT_TRUE(elemental_composition1 == elemental_composition2) << 
    "elemental_composition1 == elemental_composition2 (#1)";
  elemental_composition1["C"] = 0;
  elemental_composition2["O"] = 0;
  EXPECT_TRUE(elemental_composition1 == elemental_composition2) << 
    "elemental_composition1 == elemental_composition2 (#2)";
}

TEST_F (ElementalCompositionTestFixture, RelationalOperatorsNotEqual)
{
  zerork::ElementalComposition elemental_composition1 = zerork::ElementalComposition("AA");
  elemental_composition1["H"] = 2;
  zerork::ElementalComposition elemental_composition2 = elemental_composition1;
  elemental_composition2.name() = std::string("BB");
  EXPECT_TRUE(elemental_composition1 != elemental_composition2) << 
    "elemental_composition1 != elemental_composition2 (#1)";
  elemental_composition2 = elemental_composition1;
  elemental_composition1["H"] = 3;
  EXPECT_TRUE(elemental_composition1 != elemental_composition2) << 
    "elemental_composition1 != elemental_composition2 (#2)";
  elemental_composition2 = elemental_composition1;
  elemental_composition1["C"] = 1;
  EXPECT_TRUE(elemental_composition1 != elemental_composition2) << 
    "elemental_composition1 != elemental_composition2 (#3)";
}

TEST_F (ElementalCompositionTestFixture, RelationalOperatorsLessThan)
{
  zerork::ElementalComposition elemental_composition1 = zerork::ElementalComposition("AA");
  elemental_composition1["C"] = 2;
  elemental_composition1["H"] = 6;
  zerork::ElementalComposition elemental_composition2 = elemental_composition1;
  elemental_composition2["H"] -= 1;
  EXPECT_TRUE(elemental_composition2 < elemental_composition1) << 
    "elemental_composition1 < elemental_composition2 (#1)";
  elemental_composition2 = elemental_composition1;
  elemental_composition2.name() = "BB";
  EXPECT_TRUE(elemental_composition1 < elemental_composition2) << 
    "elemental_composition1 < elemental_composition2 (#2)";
}


TEST_F (ElementalCompositionTestFixture, Plus)
{
  zerork::ElementalComposition elemental_composition1 = zerork::ElementalComposition("AA");
  elemental_composition1["C"] = 2;
  elemental_composition1["H"] = 6;
  zerork::ElementalComposition elemental_composition2 = elemental_composition1 + elemental_composition1;
  EXPECT_EQ(elemental_composition2.GetNumTotalAtoms(), elemental_composition1.GetNumTotalAtoms()*2) << 
    "elemental_composition1 + elemental_composition1 (#1)";
  EXPECT_EQ(elemental_composition2["C"], elemental_composition1["C"]*2) << 
    "elemental_composition1 + elemental_composition1 (#2)";
  EXPECT_EQ(elemental_composition2["H"], elemental_composition1["H"]*2) << 
    "elemental_composition1 + elemental_composition1 (#3)";
}


TEST_F (ElementalCompositionTestFixture, Minus)
{
  zerork::ElementalComposition elemental_composition1 = zerork::ElementalComposition("AA");
  elemental_composition1["C"] = 2;
  elemental_composition1["H"] = 6;
  zerork::ElementalComposition elemental_composition2 = elemental_composition1 - elemental_composition1;
  EXPECT_EQ(elemental_composition2.GetNumTotalAtoms(), 0) << 
    "elemental_composition1 + elemental_composition1 (#1)";
  EXPECT_EQ(elemental_composition2["C"], 0) << 
    "elemental_composition1 + elemental_composition1 (#2)";
  EXPECT_EQ(elemental_composition2["H"], 0) << 
    "elemental_composition1 + elemental_composition1 (#3)";
}

TEST_F (ElementalCompositionTestFixture, PlusEquals)
{
  zerork::ElementalComposition elemental_composition1 = zerork::ElementalComposition("AA");
  elemental_composition1["C"] = 2;
  elemental_composition1["H"] = 6;
  zerork::ElementalComposition elemental_composition2 = elemental_composition1;
  elemental_composition2 += elemental_composition1;
  EXPECT_EQ(elemental_composition2.GetNumTotalAtoms(), elemental_composition1.GetNumTotalAtoms()*2) << 
    "elemental_composition1 + elemental_composition1 (#1)";
  EXPECT_EQ(elemental_composition2["C"], elemental_composition1["C"]*2) << 
    "elemental_composition1 + elemental_composition1 (#2)";
  EXPECT_EQ(elemental_composition2["H"], elemental_composition1["H"]*2) << 
    "elemental_composition1 + elemental_composition1 (#3)";
}


TEST_F (ElementalCompositionTestFixture, MinusEquals)
{
  zerork::ElementalComposition elemental_composition1 = zerork::ElementalComposition("AA");
  elemental_composition1["C"] = 2;
  elemental_composition1["H"] = 6;
  zerork::ElementalComposition elemental_composition2 = elemental_composition1;
  elemental_composition2 -= elemental_composition1;
  EXPECT_EQ(elemental_composition2.GetNumTotalAtoms(), 0) << 
    "elemental_composition1 + elemental_composition1 (#1)";
  EXPECT_EQ(elemental_composition2["C"], 0) << 
    "elemental_composition1 + elemental_composition1 (#2)";
  EXPECT_EQ(elemental_composition2["H"], 0) << 
    "elemental_composition1 + elemental_composition1 (#3)";
}


// --------------------------------------------------------------------------

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

