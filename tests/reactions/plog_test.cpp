#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cstdlib>

#include "zerork/mechanism.h"
#include "utilities/math_utilities.h"

const double KP1atm=1.01325e5;
const double KEnergyToTemp=4184.0/8.3144621e3;  // [K/(cal/mol)]
const double KConcentrationToKMol = 1000.0; // [mol/cm^3 -> kmol/m^3]
// test constants
const double KMinPressure    = 50.0;  // [Pa]
const double KMaxPressure    = 2.0e7; // [Pa]
const double KMinTemperature = 250.0; // [K]
const double KMaxTemperature = 3500.0; // [K]

// linear pressure dependence
double kfwd_reaction01(const double p)
{return 1.0e8*(p/KP1atm);}

// quadratic pressure dependence
double kfwd_reaction02(const double p)
{return 1024.0*(p/KP1atm)*(p/KP1atm);}

// square root pressure dependence
double kfwd_reaction03(const double p)
{return 1.6e12*sqrt(p/KP1atm);}

// inverse cubic pressure dependence
double kfwd_reaction04(const double p)
{return 8.192e-3*pow(5.0*KP1atm/p,3.0);}

// actual plog in use from tpgme
double kfwd_reaction05(const double p, const double T);

typedef struct {
  double pressure;
  double temperature;
  double mole_fraction[2];
  double concentration[2];
} TestState;

void GetRandomTestState(zerork::mechanism *mech, TestState *test);

double RunPressureTest(zerork::mechanism *mech,
                       const int fwd_rxn_id,
                       const double abs_tol,
                       const double rel_tol,
                       double (*pfunc)(double),
                       double *rel_error);

double RunPressureTemperatureTest(zerork::mechanism *mech,
                                  const int fwd_rxn_id,
                                  const double abs_tol,
                                  const double rel_tol,
                                  double (*ptfunc)(double,double),
                                  double *rel_error);

int main(int argc, char *argv[])
{
  int seed, num_samples;
  double max_abs_err,max_rel_err;
  double current_abs_err,current_rel_err;
  double test_max_abs, test_max_rel;

  if(argc != 5) {
    printf("ERROR: incorrect command line usage.\n");
    printf("       use instead %s <# sample tests/reaction>\n",
           argv[0]);
    printf("         <seed (-1 = time(0))>\n");
    printf("         <max abs. error [-]>\n");
    printf("         <max rel. error [-]>\n");
    return 1;
  }
  
  num_samples = atoi(argv[1]);
  seed        = atoi(argv[2]);
  max_abs_err = atof(argv[3]);
  max_rel_err = atof(argv[4]);

  if(seed == -1) {seed = time(0);}
  
  zerork::utilities::random01seed(seed);

  const char * TEST_MECH = "mechanisms/plog/plog_test.mech";
  const char * TEST_THERM = "mechanisms/plog/plog_test.therm";

  const char * ZERORK_DATA_DIR = std::getenv("ZERORK_DATA_DIR");
  std::string load_mech(TEST_MECH);
  std::string load_therm(TEST_THERM);
  if(ZERORK_DATA_DIR == nullptr) {
    load_mech = std::string("../../") + load_mech;
    load_therm = std::string("../../") + load_therm;
  } else {
    load_mech = std::string(ZERORK_DATA_DIR) + "/" + load_mech;
    load_therm = std::string(ZERORK_DATA_DIR) + "/" + load_therm;
  }

  zerork::mechanism mech(load_mech.c_str(), load_therm.c_str(),
                      "plog_test_parser.log");

  // reaction 1
  test_max_abs=test_max_rel=0.0;
  for(int j=0; j<num_samples; ++j) {
    current_abs_err = RunPressureTest(&mech,
                                      1, // reaction number (starting at 1)
                                      max_abs_err,
                                      max_rel_err,
                                      kfwd_reaction01,
                                      &current_rel_err);
    if(fabs(current_abs_err) > test_max_abs) {
      test_max_abs = fabs(current_abs_err);
    }
    if(fabs(current_rel_err) > test_max_rel) {
      test_max_rel = fabs(current_rel_err);
    }
  }
  printf("Reaction 1 test: rel err = %10.3e, abs err = %10.3e\n",
         test_max_rel,test_max_abs);
  // ------------------------------------------------------------------------
  // reaction 2
  test_max_abs=test_max_rel=0.0;
  for(int j=0; j<num_samples; ++j) {
    current_abs_err = RunPressureTest(&mech,
                                      2, // reaction number (starting at 1)
                                      max_abs_err,
                                      max_rel_err,
                                      kfwd_reaction02,
                                      &current_rel_err);
    if(fabs(current_abs_err) > test_max_abs) {
      test_max_abs = fabs(current_abs_err);
    }
    if(fabs(current_rel_err) > test_max_rel) {
      test_max_rel = fabs(current_rel_err);
    }
  }
  printf("Reaction 2 test: rel err = %10.3e, abs err = %10.3e\n",
         test_max_rel,test_max_abs);
  // ------------------------------------------------------------------------
  // reaction 3
  test_max_abs=test_max_rel=0.0;
  for(int j=0; j<num_samples; ++j) {
    current_abs_err = RunPressureTest(&mech,
                                      3, // reaction number (starting at 1)
                                      max_abs_err,
                                      max_rel_err,
                                      kfwd_reaction03,
                                      &current_rel_err);
    if(fabs(current_abs_err) > test_max_abs) {
      test_max_abs = fabs(current_abs_err);
    }
    if(fabs(current_rel_err) > test_max_rel) {
      test_max_rel = fabs(current_rel_err);
    }
  }
  printf("Reaction 3 test: rel err = %10.3e, abs err = %10.3e\n",
         test_max_rel,test_max_abs);
  // ------------------------------------------------------------------------
  // reaction 4
  test_max_abs=test_max_rel=0.0;
  for(int j=0; j<num_samples; ++j) {
    current_abs_err = RunPressureTest(&mech,
                                      4, // reaction number (starting at 1)
                                      max_abs_err,
                                      max_rel_err,
                                      kfwd_reaction04,
                                      &current_rel_err);
    if(fabs(current_abs_err) > test_max_abs) {
      test_max_abs = fabs(current_abs_err);
    }
    if(fabs(current_rel_err) > test_max_rel) {
      test_max_rel = fabs(current_rel_err);
    }
  }
  printf("Reaction 4 test: rel err = %10.3e, abs err = %10.3e\n",
         test_max_rel,test_max_abs);
  // reaction 5
  test_max_abs=test_max_rel=0.0;
  for(int j=0; j<num_samples; ++j) {
    current_abs_err = RunPressureTemperatureTest(&mech,
                                      5, // reaction number (starting at 1)
                                      max_abs_err,
                                      max_rel_err,
                                      kfwd_reaction05,
                                      &current_rel_err);
    if(fabs(current_abs_err) > test_max_abs) {
      test_max_abs = fabs(current_abs_err);
    }
    if(fabs(current_rel_err) > test_max_rel) {
      test_max_rel = fabs(current_rel_err);
    }
  }
  printf("Reaction 5 test: rel err = %10.3e, abs err = %10.3e\n",
         test_max_rel,test_max_abs);
  

  return 0;
}

void GetRandomTestState(zerork::mechanism *mech, TestState *test)
{
  double total_concentration;

  test->pressure = log(KMinPressure) + zerork::utilities::random01()*
    (log(KMaxPressure)-log(KMinPressure));
  test->pressure = exp(test->pressure);

  test->temperature = KMinTemperature + zerork::utilities::random01()*
    (KMaxTemperature-KMinTemperature);

  total_concentration = test->pressure/(mech->getGasConstant()*
                                        test->temperature);

  test->mole_fraction[0] = zerork::utilities::random01();
  test->mole_fraction[1] = 1.0-test->mole_fraction[0];

  test->concentration[0] = total_concentration*test->mole_fraction[0];
  test->concentration[1] = total_concentration*test->mole_fraction[1];
  
}
double RunPressureTest(zerork::mechanism *mech,
                       const int fwd_rxn_id,
                       const double abs_tol,
                       const double rel_tol,
                       double (*pfunc)(double),
                       double *rel_error)
{
  double abs_error;
  double direct_soln;
  const int num_reactions = mech->getNumReactions();
  std::vector<double> Kfwd,Krev;
  TestState random_state;

  Kfwd.assign(num_reactions,0.0);
  Krev.assign(num_reactions,0.0);

  GetRandomTestState(mech,&random_state);
  direct_soln = (*pfunc)(random_state.pressure);
  mech->getKrxnFromTC(random_state.temperature,
                      &random_state.concentration[0],
                      Kfwd.data(),
                      Krev.data());
  abs_error = Kfwd[fwd_rxn_id-1] - direct_soln;
 
  if(fabs(abs_error) <= abs_tol || fabs(abs_error/direct_soln) <= rel_tol) {
    // satisfied absolute or relative tolerance
    *rel_error = abs_error/direct_soln;
    return abs_error;
  }
  printf("Failed Both Tolerance Tests (fwd rxn id = %d):\n",fwd_rxn_id);
  printf("Absolute tolerance          : %10.4e\n",abs_tol);
  printf("Relative tolerance          : %10.4e\n",rel_tol);
  printf("  pressure           [Pa]   : %24.16e\n",random_state.pressure);
  printf("  temperature         [K]   : %24.16e\n",random_state.temperature);
  printf("  x[0] h  mole frac   [-]   : %24.16e\n",random_state.mole_fraction[0]);
  printf("  x[1] h2 mole frac   [-]   : %24.16e\n",random_state.mole_fraction[1]);
  printf("  c[0] h  conc [kmol/m^3]   : %24.16e\n",random_state.concentration[0]);
  printf("  c[1] h2 conc [kmol/m^3]   : %24.16e\n",random_state.concentration[1]);
  printf("  direct solution           : %24.16e\n",direct_soln);
  printf("  zerork::getKrxnFromTC(...): %24.16e\n", Kfwd[fwd_rxn_id-1]);
  printf("  absolute error            : %24.16e\n",abs_error);
  printf("  relative error            : %24.16e\n",abs_error/direct_soln);

  *rel_error = abs_error/direct_soln;
  return abs_error;
}

// actual plog in use from tpgme mechanism
double kfwd_reaction05(const double p, const double T)
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
    id = 0;
  }
  else if(p_atm >= p_pts[n_pts-1]) {
    id = n_pts -2;
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
double RunPressureTemperatureTest(zerork::mechanism *mech,
                                  const int fwd_rxn_id,
                                  const double abs_tol,
                                  const double rel_tol,
                                  double (*ptfunc)(double,double),
                                  double *rel_error)
{
  double abs_error;
  double direct_soln;
  const int num_reactions = mech->getNumReactions();
  std::vector<double> Kfwd,Krev;
  TestState random_state;

  Kfwd.assign(num_reactions,0.0);
  Krev.assign(num_reactions,0.0);

  GetRandomTestState(mech,&random_state);
  direct_soln = (*ptfunc)(random_state.pressure,
                          random_state.temperature);
  mech->getKrxnFromTC(random_state.temperature,
                      &random_state.concentration[0],
                      Kfwd.data(),
                      Krev.data());
  abs_error = Kfwd[fwd_rxn_id-1] - direct_soln;
 
  if(fabs(abs_error) <= abs_tol || fabs(abs_error/direct_soln) <= rel_tol) {
    // satisfied absolute or relative tolerance
    *rel_error = abs_error/direct_soln;
    return abs_error;
  }
  printf("Failed Both Tolerance Tests (fwd rxn id = %d):\n",fwd_rxn_id);
  printf("Absolute tolerance          : %10.4e\n",abs_tol);
  printf("Relative tolerance          : %10.4e\n",rel_tol);
  printf("  pressure           [Pa]   : %24.16e\n",random_state.pressure);
  printf("  temperature         [K]   : %24.16e\n",random_state.temperature);
  printf("  x[0] h  mole frac   [-]   : %24.16e\n",random_state.mole_fraction[0]);
  printf("  x[1] h2 mole frac   [-]   : %24.16e\n",random_state.mole_fraction[1]);
  printf("  c[0] h  conc [kmol/m^3]   : %24.16e\n",random_state.concentration[0]);
  printf("  c[1] h2 conc [kmol/m^3]   : %24.16e\n",random_state.concentration[1]);
  printf("  direct solution           : %24.16e\n",direct_soln);
  printf("  zerork::getKrxnFromTC(...): %24.16e\n", Kfwd[fwd_rxn_id-1]);
  printf("  absolute error            : %24.16e\n",abs_error);
  printf("  relative error            : %24.16e\n",abs_error/direct_soln);

  *rel_error = abs_error/direct_soln;
  return abs_error;
}
