#include <math.h>
#include <string.h> // memcpy

#include <vector>
#include <iostream> // for file reading
#include <fstream>  // for file reading

#include <gtest/gtest.h>

#include <math_utilities.h>
#include <string_utilities.h>
#include <file_utilities.h>

const int NUM_SIMPLE = 11;
const int SIMPLE_I[] = {1024, 1, 512, 2, 256, 4, 128, 8, 64, 16, 32};
const double SIMPLE_D[] = {-1024.0, 1.0, -512.0,  2.0, -256.0, 4.0, -128.0,
		           8.0, -64.0, 16.0, -32.0};
const int NUM_DUPLICATE = 14;
const double DUPLICATE_D[] = {-1024.0, 1.0, -512.0,  2.0, -256.0, 4.0, -128.0,
		              8.0, -64.0, 16.0, -32.0, -512.0, 2.0, -64.0};
const int    VALUES_I[] = {0,1,2,3,4,5,6,7,8,9,10};
const double VALUES_D[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};

const int NUM_INTERP = 8;
const double X_INTERP[] = {-1.5, -1.0, -0.5, 0.0, 1.0, 2.0, 3.0, 4.0};
const double F_INTERP[] = { 3.0,  2.0,  1.0, 0.0, 1.0, 4.0, 9.0,16.0};

const double X_SCRAMBLE[] = { -1.0, 3.0,  4.0, -0.5, 1.0, -1.5, 0.0, 2.0};
const double F_SCRAMBLE[] = {  2.0, 9.0, 16.0,  1.0, 1.0,  3.0, 0.0, 4.0};
const double X_DUPLICATE[] = { 1.0, 2.0, 3.0, -1.5, -1.0, -0.5, 0.0, 1.0, -1.0, 2.0, 3.0, -0.5, 0.0,  4.0};
const double F_DUPLICATE[] = { 1.0, 4.0, 9.0,  3.0,  2.0,  1.0, 0.0, 1.0,  2.0, 4.0, 9.0,  1.0, 0.0, 16.0};

// TEST (SortVectors, SimpleInt)
// {
//   std::vector<int> keys, values;
//   keys.assign(NUM_SIMPLE,0);
//   values.assign(NUM_SIMPLE,0);

//   memcpy(&keys[0],   SIMPLE_I, NUM_SIMPLE*sizeof(int));
//   memcpy(&values[0], VALUES_I, NUM_SIMPLE*sizeof(int));

//   zerork::utilities::SortVectors(keys, values);

//   EXPECT_EQ(keys[ 0],    1);
//   EXPECT_EQ(keys[ 1],    2);
//   EXPECT_EQ(keys[ 2],    4);
//   EXPECT_EQ(keys[ 3],    8);
//   EXPECT_EQ(keys[ 4],   16);
//   EXPECT_EQ(keys[ 5],   32);
//   EXPECT_EQ(keys[ 6],   64);
//   EXPECT_EQ(keys[ 7],  128);
//   EXPECT_EQ(keys[ 8],  256);
//   EXPECT_EQ(keys[ 9],  512);
//   EXPECT_EQ(keys[10], 1024);
// }

TEST (SortVectors, SimpleDouble)
{
  std::vector<double> keys, values;

  keys.assign(NUM_SIMPLE,0);
  values.assign(NUM_SIMPLE,0);
  memcpy(&keys[0],   SIMPLE_D, NUM_SIMPLE*sizeof(double));
  memcpy(&values[0], VALUES_D, NUM_SIMPLE*sizeof(double));

  zerork::utilities::SortVectors(&keys, &values);

  EXPECT_EQ(keys[ 0],-1024.);
  EXPECT_EQ(keys[ 1], -512.);
  EXPECT_EQ(keys[ 2], -256.);
  EXPECT_EQ(keys[ 3], -128.);
  EXPECT_EQ(keys[ 4],  -64.);
  EXPECT_EQ(keys[ 5],  -32.);
  EXPECT_EQ(keys[ 6],    1.);
  EXPECT_EQ(keys[ 7],    2.);
  EXPECT_EQ(keys[ 8],    4.);
  EXPECT_EQ(keys[ 9],    8.);
  EXPECT_EQ(keys[10],   16.);
  EXPECT_EQ(values[ 0], 0.);
  EXPECT_EQ(values[ 1], 2.);
  EXPECT_EQ(values[ 2], 4.);
  EXPECT_EQ(values[ 3], 6.);
  EXPECT_EQ(values[ 4], 8.);
  EXPECT_EQ(values[ 5],10.);
  EXPECT_EQ(values[ 6], 1.);
  EXPECT_EQ(values[ 7], 3.);
  EXPECT_EQ(values[ 8], 5.);
  EXPECT_EQ(values[ 9], 7.);
  EXPECT_EQ(values[10], 9.);

}
TEST (SortVectors, DuplicateDouble)
{
  std::vector<double> keys, values;

  keys.assign(NUM_DUPLICATE,0);
  values.assign(NUM_DUPLICATE,0);
  memcpy(&keys[0],   DUPLICATE_D, NUM_DUPLICATE*sizeof(double));
  memcpy(&values[0], VALUES_D, NUM_DUPLICATE*sizeof(double));

  zerork::utilities::SortVectors(&keys, &values);
  EXPECT_EQ(keys[ 0],-1024.); EXPECT_EQ(values[ 0], 0.);
  EXPECT_EQ(keys[ 1], -512.);
  EXPECT_EQ(keys[ 2], -512.);
  EXPECT_EQ(keys[ 3], -256.); EXPECT_EQ(values[ 3], 4.);
  EXPECT_EQ(keys[ 4], -128.); EXPECT_EQ(values[ 4], 6.);
  EXPECT_EQ(keys[ 5],  -64.);
  EXPECT_EQ(keys[ 6],  -64.);
  EXPECT_EQ(keys[ 7],  -32.); EXPECT_EQ(values[ 7], 10.);
  EXPECT_EQ(keys[ 8],    1.); EXPECT_EQ(values[ 8], 1.);
  EXPECT_EQ(keys[ 9],    2.);
  EXPECT_EQ(keys[10],    2.);
  EXPECT_EQ(keys[11],    4.); EXPECT_EQ(values[11], 5.);
  EXPECT_EQ(keys[12],    8.); EXPECT_EQ(values[12], 7.);
  EXPECT_EQ(keys[13],   16.); EXPECT_EQ(values[13], 9.);

  // Duplicates
  EXPECT_TRUE(values[ 1] == 2. || values[ 1] == 11.) <<
    "values[ 1] = " << values[ 1] <<std:: endl;
  EXPECT_TRUE(values[ 2] == 2. || values[ 2] == 11.) <<
    "values[ 1] = " << values[ 2] << std::endl;
  EXPECT_TRUE(values[ 5] == 8. || values[ 5] == 13.) <<
    "values[ 1] = " << values[ 5] << std::endl;
  EXPECT_TRUE(values[ 6] == 8. || values[ 6] == 13.) <<
    "values[ 1] = " << values[ 6] << std::endl;
  EXPECT_TRUE(values[ 9] == 3. || values[ 9] == 12.) <<
    "values[ 1] = " << values[ 9] << std::endl;
  EXPECT_TRUE(values[10] == 3. || values[10] == 12.) <<
    "values[ 1] = " << values[10] << std::endl;
}


TEST(InterpolationTable, OrderedListNoExtrapolation)
{
  zerork::utilities::InterpolationTable *table=NULL;
  table = new zerork::utilities::InterpolationTable(NUM_INTERP,
                                            X_INTERP,
                                            F_INTERP,
                                            zerork::utilities::LINEAR,
                                            false);
  ASSERT_TRUE(table != NULL) << "Constructor failed" << std::endl;
  EXPECT_DOUBLE_EQ(table->Interpolate(-10.0),   3.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.5),    3.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.25),   2.5);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.9),    1.8);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.5),    1.0);
  EXPECT_NEAR(table->Interpolate(-1.0e-4), 2.0e-4, 1.0e-16);
  EXPECT_DOUBLE_EQ(table->Interpolate(0.3),     0.3);
  EXPECT_DOUBLE_EQ(table->Interpolate(1.7),     3.1);
  EXPECT_DOUBLE_EQ(table->Interpolate(2.2),     5.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.0),     9.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.4),    11.8);
  EXPECT_DOUBLE_EQ(table->Interpolate(4.0),    16.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(10.0),   16.0);

  delete table;
}
TEST(InterpolationTable, OrderedListExtrapolation)
{
  zerork::utilities::InterpolationTable *table=NULL;
  table = new zerork::utilities::InterpolationTable(NUM_INTERP,
                                            X_INTERP,
                                            F_INTERP,
                                            zerork::utilities::LINEAR,
                                            true);
  ASSERT_TRUE(table != NULL) << "Constructor failed" << std::endl;
  EXPECT_DOUBLE_EQ(table->Interpolate(-10.0),  20.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.5),    3.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.25),   2.5);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.9),    1.8);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.5),    1.0);
  EXPECT_NEAR(table->Interpolate(-1.0e-4), 2.0e-4, 1.0e-16);
  EXPECT_DOUBLE_EQ(table->Interpolate(0.3),     0.3);
  EXPECT_DOUBLE_EQ(table->Interpolate(1.7),     3.1);
  EXPECT_DOUBLE_EQ(table->Interpolate(2.2),     5.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.0),     9.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.4),    11.8);
  EXPECT_DOUBLE_EQ(table->Interpolate(4.0),    16.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(10.0),   58.0);

  delete table;
}


TEST(InterpolationTable, UnsortedListNoExtrapolation)
{
  zerork::utilities::InterpolationTable *table=NULL;
  table = new zerork::utilities::InterpolationTable(NUM_INTERP,
                                            X_SCRAMBLE,
                                            F_SCRAMBLE,
                                            zerork::utilities::LINEAR,
                                            false);
  ASSERT_TRUE(table != NULL) << "Constructor failed" << std::endl;
  EXPECT_DOUBLE_EQ(table->Interpolate(-10.0),   3.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.5),    3.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.25),   2.5);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.9),    1.8);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.5),    1.0);
  EXPECT_NEAR(table->Interpolate(-1.0e-4), 2.0e-4, 1.0e-16);
  EXPECT_DOUBLE_EQ(table->Interpolate(0.3),     0.3);
  EXPECT_DOUBLE_EQ(table->Interpolate(1.7),     3.1);
  EXPECT_DOUBLE_EQ(table->Interpolate(2.2),     5.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.0),     9.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.4),    11.8);
  EXPECT_DOUBLE_EQ(table->Interpolate(4.0),    16.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(10.0),   16.0);

  delete table;
}

TEST(InterpolationTable, UnsortedListExtrapolation)
{
  zerork::utilities::InterpolationTable *table=NULL;
  table = new zerork::utilities::InterpolationTable(NUM_INTERP,
                                            X_SCRAMBLE,
                                            F_SCRAMBLE,
                                            zerork::utilities::LINEAR,
                                            true);
  ASSERT_TRUE(table != NULL) << "Constructor failed" << std::endl;
  EXPECT_DOUBLE_EQ(table->Interpolate(-10.0),  20.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.5),    3.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.25),   2.5);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.9),    1.8);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.5),    1.0);
  EXPECT_NEAR(table->Interpolate(-1.0e-4), 2.0e-4, 1.0e-16);
  EXPECT_DOUBLE_EQ(table->Interpolate(0.3),     0.3);
  EXPECT_DOUBLE_EQ(table->Interpolate(1.7),     3.1);
  EXPECT_DOUBLE_EQ(table->Interpolate(2.2),     5.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.0),     9.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.4),    11.8);
  EXPECT_DOUBLE_EQ(table->Interpolate(4.0),    16.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(10.0),   58.0);

  delete table;
}
TEST(InterpolationTable, DuplicateListNoExtrapolation)
{
  zerork::utilities::InterpolationTable *table=NULL;
  table = new zerork::utilities::InterpolationTable(NUM_DUPLICATE,
                                            X_DUPLICATE,
                                            F_DUPLICATE,
                                            zerork::utilities::LINEAR,
                                            false);
  ASSERT_TRUE(table != NULL) << "Constructor failed" << std::endl;
  EXPECT_DOUBLE_EQ(table->Interpolate(-10.0),   3.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.5),    3.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.25),   2.5);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.9),    1.8);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.5),    1.0);
  EXPECT_NEAR(table->Interpolate(-1.0e-4), 2.0e-4, 1.0e-16);
  EXPECT_DOUBLE_EQ(table->Interpolate(0.3),     0.3);
  EXPECT_DOUBLE_EQ(table->Interpolate(1.7),     3.1);
  EXPECT_DOUBLE_EQ(table->Interpolate(2.2),     5.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.0),     9.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.4),    11.8);
  EXPECT_DOUBLE_EQ(table->Interpolate(4.0),    16.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(10.0),   16.0);

  delete table;
}

TEST(InterpolationTable, DuplicateListExtrapolation)
{
  zerork::utilities::InterpolationTable *table=NULL;
  table = new zerork::utilities::InterpolationTable(NUM_DUPLICATE,
                                            X_DUPLICATE,
                                            F_DUPLICATE,
                                            zerork::utilities::LINEAR,
                                            true);
  ASSERT_TRUE(table != NULL) << "Constructor failed" << std::endl;
  EXPECT_DOUBLE_EQ(table->Interpolate(-10.0),  20.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.5),    3.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.25),   2.5);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.9),    1.8);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.5),    1.0);
  EXPECT_NEAR(table->Interpolate(-1.0e-4), 2.0e-4, 1.0e-16);
  EXPECT_DOUBLE_EQ(table->Interpolate(0.3),     0.3);
  EXPECT_DOUBLE_EQ(table->Interpolate(1.7),     3.1);
  EXPECT_DOUBLE_EQ(table->Interpolate(2.2),     5.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.0),     9.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.4),    11.8);
  EXPECT_DOUBLE_EQ(table->Interpolate(4.0),    16.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(10.0),   58.0);

  delete table;
}

TEST(InterpolationTable, OrderedListNoExtrapolationCubic)
{
  //N.B. Interpolated "EXPECT_NEAR" values generated from
  // spot check with comparison to GSL, not numerical analysis
  zerork::utilities::InterpolationTable *table=NULL;
  table = new zerork::utilities::InterpolationTable(NUM_INTERP,
                                            X_INTERP,
                                            F_INTERP,
                                            zerork::utilities::CUBIC_SPLINE,
                                            false);
  ASSERT_TRUE(table != NULL) << "Constructor failed" << std::endl;
  EXPECT_DOUBLE_EQ(table->Interpolate(-10.0),   3.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.5),    3.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(-0.5),    1.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(0.0),     0.0);
  EXPECT_NEAR(table->Interpolate(0.5),  0.0676248, 1.0e-5);
  EXPECT_DOUBLE_EQ(table->Interpolate(2.0),     4.0);
  EXPECT_NEAR(table->Interpolate(2.5),  6.21227, 1.0e-5);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.0),     9.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(4.0),    16.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(10.0),   16.0);

  delete table;
}

TEST(InterpolationTable, OrderedListExtrapolationCubic)
{
  //N.B. Interpolated "EXPECT_NEAR" values generated from
  // spot check with comparison to GSL, not numerical analysis
  zerork::utilities::InterpolationTable *table=NULL;
  table = new zerork::utilities::InterpolationTable(NUM_INTERP,
                                            X_INTERP,
                                            F_INTERP,
                                            zerork::utilities::CUBIC_SPLINE,
                                            true);
  ASSERT_TRUE(table != NULL) << "Constructor failed" << std::endl;
  EXPECT_NEAR(table->Interpolate(-5.0),  4.459, 1.0e-4);
  EXPECT_DOUBLE_EQ(table->Interpolate(-1.5),    3.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(0.0),     0.0);
  EXPECT_NEAR(table->Interpolate(0.5),  0.0676248, 1.0e-5);
  EXPECT_DOUBLE_EQ(table->Interpolate(2.0),     4.0);
  EXPECT_NEAR(table->Interpolate(2.5),  6.21227, 1.0e-5);
  EXPECT_DOUBLE_EQ(table->Interpolate(3.0),     9.0);
  EXPECT_DOUBLE_EQ(table->Interpolate(4.0),    16.0);
  EXPECT_NEAR(table->Interpolate(5.0),  23, 1.0e-4);

  delete table;
}

TEST(InterpolationTable, OrderedListNoExtrapolationFirstDerivative)
{
  zerork::utilities::InterpolationTable *table=NULL;
  table = new zerork::utilities::InterpolationTable(NUM_INTERP,
                                            X_INTERP,
                                            F_INTERP,
                                            zerork::utilities::LINEAR,
                                            false);
  ASSERT_TRUE(table != NULL) << "Constructor failed" << std::endl;
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(-10.0),   0.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(-1.25),   -2.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(-0.75),   -2.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(-0.25),   -2.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(0.5),   1.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(1.5),   3.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(2.5),   5.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(3.5),   7.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(10.0),   0.0);

  delete table;
}

TEST(InterpolationTable, OrderedListExtrapolationFirstDerivative)
{
  zerork::utilities::InterpolationTable *table=NULL;
  table = new zerork::utilities::InterpolationTable(NUM_INTERP,
                                            X_INTERP,
                                            F_INTERP,
                                            zerork::utilities::LINEAR,
                                            true);
  ASSERT_TRUE(table != NULL) << "Constructor failed" << std::endl;
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(-10.0),   -2.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(-1.25),   -2.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(-0.75),   -2.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(-0.25),   -2.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(0.5),   1.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(1.5),   3.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(2.5),   5.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(3.5),   7.0);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(10.0),   7.0);

  delete table;
}

TEST(InterpolationTable, OrderedListNoExtrapolationCubicFirstDerivative)
{
  zerork::utilities::InterpolationTable *table=NULL;
  table = new zerork::utilities::InterpolationTable(NUM_INTERP,
                                            X_INTERP,
                                            F_INTERP,
                                            zerork::utilities::CUBIC_SPLINE,
                                            false);
  ASSERT_TRUE(table != NULL) << "Constructor failed" << std::endl;
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(-5.0),   0.0);
  EXPECT_NEAR(table->InterpolateFirstDerivative(-1.2),   -1.99736, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative(-0.7),   -1.99340, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative(-0.2),   -2.02902, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative( 0.5),    1.20648, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative( 1.5),    2.95179, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative( 2.5),    4.98637, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative( 3.5),    7.10273, 1.0e-5);
  EXPECT_DOUBLE_EQ(table->InterpolateFirstDerivative(5.0),   0.0);

  delete table;
}

TEST(InterpolationTable, OrderedListExtrapolationCubicFirstDerivative)
{
  //N.B. Interpolated "EXPECT_NEAR" values generated from
  // spot check with comparison to GSL, not numerical analysis
  zerork::utilities::InterpolationTable *table=NULL;
  table = new zerork::utilities::InterpolationTable(NUM_INTERP,
                                            X_INTERP,
                                            F_INTERP,
                                            zerork::utilities::CUBIC_SPLINE,
                                            true);
  ASSERT_TRUE(table != NULL) << "Constructor failed" << std::endl;
  EXPECT_NEAR(table->InterpolateFirstDerivative(-5.0),    2.81539, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative(-1.2),   -1.99736, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative(-0.7),   -1.99340, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative(-0.2),   -2.02902, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative( 0.5),    1.20648, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative( 1.5),    2.95179, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative( 2.5),    4.98637, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative( 3.5),    7.10273, 1.0e-5);
  EXPECT_NEAR(table->InterpolateFirstDerivative( 5.0),    6.17820, 1.0e-5);

  delete table;
}

TEST(erfc_inv, AccuracyCheck)
{
  std::ifstream input_file("data/erfc_inv_test.dat");
  std::string line;
  int line_num = 0;
  std::vector<std::string> token_vector;
  while(zerork::utilities::GetAnyLine(input_file,&line)) {
    size_t num_tokens =
      zerork::utilities::SplitStringToVector(line, zerork::utilities::WHITESPACE, &token_vector);
    if(num_tokens != 2) continue; //EOF

    double q = std::stod(token_vector[0]);
    double x = std::stod(token_vector[1]);
    double calc_x = zerork::utilities::erfc_inv(q);
    EXPECT_NEAR(x, calc_x, 5.0e-3);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
