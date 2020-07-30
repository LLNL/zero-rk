#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "distribution.h"

void fill_quadratic(const int j_start,
                    const int j_end,
                    zerork::utilities::Distribution *dist);
void fill_linear(const int j_start,
                 const int j_end,
                 zerork::utilities::Distribution *dist);

const int num_samples_test01 = 41;
const int num_bins_test01 = 10;
const double min_test01 = 1.0;
const double max_test01 = 1024.0;
bool checkSuccessTest01(zerork::utilities::Distribution *dist);

const int num_bins_test02 = 17;
const double test02_soln = 13.0;
const double min_test02 = 0.5;
const double max_test02 = min_test02+num_bins_test02*test02_soln;
bool checkSuccessTest02(zerork::utilities::Distribution *dist);

int main(int argc, char *argv[])
{
  zerork::utilities::Distribution dist01(num_bins_test01,
                           min_test01,
                           max_test01,
                           true);
  zerork::utilities::Distribution *dist01a;
  fill_quadratic(0,num_samples_test01-1,&dist01);
  dist01a = new zerork::utilities::Distribution(dist01); // new copy constructor
  zerork::utilities::Distribution dist01b(dist01);  // copy constructor
  zerork::utilities::Distribution dist01c; // default constructor (unit interval, 1 bin)

  printf("Test 1 : quadratic integers(+0.5) in log-spaced bins\n"); 
  if(checkSuccessTest01(&dist01)) {
    printf("  SUCCESS\n");
  }
  printf("Test 1a: use new Distribution(&Distribution) copy constructor\n"); 
  if(checkSuccessTest01(dist01a)) {
    printf("  SUCCESS\n");
  }  
  printf("Test 1b: use direct Distribution(&Distribution) copy constructor\n"); 
  if(checkSuccessTest01(&dist01b)) {
    printf("  SUCCESS\n");
  }
  dist01c = dist01; // use operator=
  printf("Test 1c: use default constructor with assignment (=) operator\n"); 
  if(checkSuccessTest01(&dist01c)) {
    printf("  SUCCESS\n");
  }

  // perform a two partial fills
  dist01a->zeroBinWeights();
  dist01b.zeroBinWeights();
  fill_quadratic(0,10,dist01a);
  fill_quadratic(11,num_samples_test01-1,&dist01b);
  printf("Test 1d: use zeroBinWeights, assign 2 partial distributions,\n");
  printf("         and then use the += operator for the full distribution.\n");
  dist01b+=(*dist01a);
  if(checkSuccessTest01(&dist01b)) {
    printf("  SUCCESS\n");
  }
  // perform a three partial fills
  dist01a->zeroBinWeights();
  dist01b.zeroBinWeights();
  dist01c.zeroBinWeights();
  fill_quadratic(0,7,dist01a);
  fill_quadratic(8,17,&dist01b);
  fill_quadratic(18,num_samples_test01-1,&dist01c);
  printf("Test 1e: use zeroBinWeights, assign 3 partial distributions,\n");
  printf("         and then use ( = a+b+c) operators for the full distribution.\n");
  dist01 = (*dist01a) + dist01b + dist01c;
  if(checkSuccessTest01(&dist01)) {
    printf("  SUCCESS\n");
  }
  

  zerork::utilities::Distribution dist02(num_bins_test02,
                           min_test02,
                           max_test02,
                           false);
  zerork::utilities::Distribution *dist02a, dist02c;
  fill_linear(-10,250,&dist02);

  dist02a = new zerork::utilities::Distribution(dist02); // new copy constructor
  zerork::utilities::Distribution dist02b(dist02);  // copy constructor

  printf("Test 2 : linear integers in linear-spaced bins\n"); 
  if(checkSuccessTest02(&dist02)) {
    printf("  SUCCESS\n");
  }
  printf("Test 2a: use new Distribution(&Distribution) copy constructor\n"); 
  if(checkSuccessTest02(dist02a)) {
    printf("  SUCCESS\n");
  }  
  printf("Test 2b: use direct Distribution(&Distribution) copy constructor\n"); 
  if(checkSuccessTest02(&dist02b)) {
    printf("  SUCCESS\n");
  }  
  dist02c = dist02; // use operator=
  printf("Test 2c: use default constructor with assignment (=) operator\n"); 
  if(checkSuccessTest02(&dist02c)) {
    printf("  SUCCESS\n");
  }
  // perform a two partial fills
  dist02a->zeroBinWeights();
  dist02b.zeroBinWeights();
  fill_linear(-10,77,dist02a);
  fill_linear(78,250,&dist02b);
  printf("Test 2d: use zeroBinWeights, assign 2 partial distributions,\n");
  printf("         and then use the += operator for the full distribution.\n");
  dist02b+=(*dist02a);
  if(checkSuccessTest02(&dist02b)) {
    printf("  SUCCESS\n");
  }
  // perform a three partial fills
  dist02a->zeroBinWeights();
  dist02b.zeroBinWeights();
  dist02c.zeroBinWeights();
  fill_linear(-10,77,dist02a);
  fill_linear(111,250,&dist02b);
  fill_linear(78,110,&dist02c);
  printf("Test 2e: use zeroBinWeights, assign 3 partial distributions,\n");
  printf("         and then use ( = a+b+c) operators for the full distribution.\n");
  dist02 = (*dist02a) + dist02b + dist02c;
  if(checkSuccessTest02(&dist02)) {
    printf("  SUCCESS\n");
  }

  printf("Test 3a: use the += on two different distributions\n");
  dist01b+=dist02b;
    

  

  delete dist01a;
  delete dist02a;

  return 0;
}

void fill_quadratic(const int j_start,
                    const int j_end,
                    zerork::utilities::Distribution *dist)
{
  double val;
  for(int j=j_start; j<= j_end; ++j) {
    val=static_cast<double>(j);
    val*=val;
    val+=0.5;
    dist->addValue(val,1.0);
    //int ret_id = dist->addValue(val,1.0);
    //printf("val=%.18g, ret id=%d, bin id=%d\n",
    //       val,ret_id,dist->getBinIndex(val));
  }
}
void fill_linear(const int j_start,
                 const int j_end,
                 zerork::utilities::Distribution *dist)
{
  double val;
  for(int j=j_start; j<= j_end; ++j) {
    val=static_cast<double>(j);
    dist->addValue(val,1.0);
  }
}
// Test01 checks for a log distribution in power-of-2 bins from 1 to 1024,
// for 41 samples j=0,1,...,40 of the function f(j) = j*j+0.5.
bool checkSuccessTest01(zerork::utilities::Distribution *dist)
{
  bool found_fail=false;
  const double under_test01 = 1; // j < 1 
  const double over_test01 = 9; // j >= 32
  int test01_soln[num_bins_test01];
  test01_soln[0] = 1; // 1.5 quadratic integers+0.5 in [1,2)
  test01_soln[1] = 0; // --- quadratic integers+0.5 in [2,4)
  test01_soln[2] = 1; // 4.5 quadratic integers+0.5 in [4,8)
  test01_soln[3] = 1; // 9.5 quadratic integers+0.5 in [8,16)
  test01_soln[4] = 2; // 16.5, 25.5 quadratic integers+0.5 in [16,32)
  test01_soln[5] = 2; // 36.5, 49.5 quadratic integers+0.5 in [32,64)
  test01_soln[6] = 4; // 64.5, 81.5, 100.5, 121.5 in [64,128)
  test01_soln[7] = 4; // 144.5, 169.5, 196.5, 225.5 in [128,256)
  test01_soln[8] = 7; // 16, 17, 18, 19, 20, 21, 22**2+0.5 = 484.5
  test01_soln[9] = 9; // 23, ..., 31**2+0.5 = 961.5

  if(dist->num_bins() != num_bins_test01) {
    printf("**FAILED - no. distribution bins %d != %d\n",
	   dist->num_bins(),
           num_bins_test01);
    return false;
  }

  for(int j=0; j<num_bins_test01; j++) {
    if(fabs(static_cast<double>(test01_soln[j])-
            dist->getBinWeight(j)) > 0.1) {
  
      printf("**FAILED - "); 
      printf("Bin [%2d]: in [%4g, %4g) distribution count %4g - exact %d\n",
             j,
             dist->getBinMin(j),
             dist->getBinMax(j),
             dist->getBinWeight(j),
             test01_soln[j]);
      found_fail = true;
    }
  }
  if(fabs(under_test01-dist->under_range_weight()) > 0.1) {
    printf("**FAILED - ");
    printf("Under range: (-\\infty, %4g) distribution count %4g - exact %4g\n",
           dist->min_range(),
           dist->under_range_weight(),
           under_test01);
    found_fail = true;
  }
  if(fabs(over_test01-dist->over_range_weight()) > 0.1) {
    printf("**FAILED - ");
    printf("Over  range: [%4g, \\infty) distribution count %4g - exact %4g\n",
           dist->max_range(),
           dist->over_range_weight(),
           over_test01);
    found_fail = true;
  }
  // return if test was successful
  if(found_fail) {
    return false;
  } else {
    return true;
  }
}
// Test02 checks for a linear distribution in 17 equally spaced bins from
// 0.5 to 17*13+0.5 for 261 samples j=-10,-9,...,250 of the function
// f(j) = j.
bool checkSuccessTest02(zerork::utilities::Distribution *dist)
{
  bool found_fail=false;
  const double under_test02 = 11; // -10,-9,...,0 
  const double over_test02 = 29;  // 222,...,250  

  if(dist->num_bins() != num_bins_test02) {
    printf("**FAILED - no. distribution bins %d != %d\n",
	   dist->num_bins(),
           num_bins_test02);
    return false;
  }

  for(int j=0; j<num_bins_test02; ++j) {
    if(fabs(test02_soln-dist->getBinWeight(j)) > 0.1) {
      printf("**FAILED - "); 
      printf("Bin [%2d]: in [%4g, %4g) distribution count %4g - exact %g\n",
             j,
             dist->getBinMin(j),
             dist->getBinMax(j),
             dist->getBinWeight(j),
             test02_soln);
      found_fail = true;
    }
  }
  if(fabs(under_test02-dist->under_range_weight()) > 0.1) {
    printf("**FAILED - ");
    printf("Under range: (-\\infty, %4g) distribution count %4g - exact %4g\n",
           dist->min_range(),
           dist->under_range_weight(),
           under_test02);
    found_fail = true;
  }
  if(fabs(over_test02-dist->over_range_weight()) > 0.1) {
    printf("**FAILED - ");
    printf("Over  range: [%4g, \\infty) distribution count %4g - exact %4g\n",
           dist->max_range(),
           dist->over_range_weight(),
           over_test02);
    found_fail = true;
  }
  // return if test was successful
  if(found_fail) {
    return false;
  } else {
    return true;
  }
}
