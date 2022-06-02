#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <string>
#include <vector>

#include <zerork/constants.h>
#include <reactor/const_pressure_reactor.h>
#include <utilities/math_utilities.h>
#include <utilities/file_utilities.h>

using zerork::utilities::null_filename;

// TEST CONSTANTS
const double MIN_PRESSURE    = 1.0e5;  // [Pa]
const double MAX_PRESSURE    = 1.0e7;  // [Pa]

const double MIN_TEMPERATURE =  300.0; // [K]
const double MAX_TEMPERATURE = 3500.0; // [K]

const double REF_TEMPERATURE = 3000.0; // [K]

const double TEST_RTOL = 1.0e-10;
const double TEST_ATOL = 1.0e-20;

// Mechanism data constants
const char TEST_MECH[]  = "mechanisms/ideal/hydrogen_recombination_rev_zero.mech";
const char TEST_THERM[] = "mechanisms/ideal/const_specific_heat.therm";
const int NUM_SPECIES = 2;
const int NUM_STATES  = 4;
const double RATE_CONST = 2.0e9; // needs to be 1000 times smaller than 
                                 // chemkin mecahnism value to convert from
                                 // the rate constant units:
                                 // 1 [cm^3/mol/s] = 1000 [cm^3/kmol/s]
                                 //                = 0.001 [m^3/kmol/s]

const double ATOMIC_MASS =  1.00794;   // [kg/kmol] note value differs with
                                       // cantera 1.7
const double ATOMIC_A0   =  2.5;
const double ATOMIC_A5   =  2.54737667e+04;
const double DIATOMIC_A0 =  3.5;
const double DIATOMIC_A5 = -1.04352500e+03; 

static int GetExactDerivative(const double atomic_mass_fraction,
                              const double temperature,
                              const double pressure,
                              std::vector<double> *derivative);

static int GetReactorDerivative(const double atomic_mass_fraction,
                                const double temperature,
                                const double pressure,
                                std::vector<double> *state,
                                std::vector<double> *derivative);

static int RandomInt(const int a, const int b);

static void RandomNormalizedVector(const size_t num_elements,
                                   const size_t num_digits, 
                                   std::vector<double> *v);

bool SameVector(const std::vector<double>A, 
                const std::vector<double>B,
                const double rtol,
                const double atol);


int main(int argc, char *argv[])
{
  int seed = -1;     // default indicates time(0) used as seed
  int num_runs = 1;
  int num_passed = 0;
  std::string outcome;
  double atomic_mass_fraction, temperature;
  double pressure=1.0e5;
  std::vector<double> reactor_state;
  std::vector<double> reactor_derivative;
  std::vector<double> exact_derivative;

  if(argc == 4) {
    // command line contains the atomic mass fraction, temperature and
    // pressure
    atomic_mass_fraction = atof(argv[1]);
    temperature          = atof(argv[2]);
    pressure             = atof(argv[3]);
  }

  if(argc == 3) {
    // command line contains num_runs and seed
    seed = atoi(argv[2]);
  } 
  if(argc == 2 || argc == 3) {
    num_runs = atoi(argv[1]);
  }
  if(num_runs <= 0) {
    num_runs = 1;
  }
  if(seed == -1) {
    seed = time(0);
  }
  zerork::utilities::random01seed(seed);

  for(int j=0; j<num_runs; ++j) {
    if(argc != 4) {
      temperature = (double)RandomInt((int)MIN_TEMPERATURE,
                                      (int)MAX_TEMPERATURE);
      pressure = (double)RandomInt((int)MIN_PRESSURE,
                                   (int)MAX_PRESSURE);
      atomic_mass_fraction = 1.0e-4*RandomInt(0,10000);    
    }
    GetReactorDerivative(atomic_mass_fraction,
                         temperature,
                         pressure,
                         &reactor_state,
                         &reactor_derivative);
    GetExactDerivative(atomic_mass_fraction,
                       temperature,
                       pressure,
                       &exact_derivative);

    if(SameVector(exact_derivative,reactor_derivative,TEST_RTOL,TEST_ATOL)) {
      outcome = "PASSED";
      ++num_passed;
    } else {
      outcome = "FAILED";
    }
    printf("[%s] RUN %2d: %s\n",
           argv[0],
           j+1,
           outcome.c_str());
    if(outcome == "FAILED" || argc == 4) {
      printf("    Initial State:\n");
      printf("        T    = %.10g [K]\n",temperature);
      printf("        p    = %.10g [Pa]\n",pressure);
      printf("        y[H] = %.10g [-]\n",atomic_mass_fraction);
      for(int k=0; k<NUM_STATES; ++k) {
        printf("    y[%2d]: %.10g  dy/dt(reactor) = %16.9e  dy/dt(exact) = %16.9e\n",
               k+1,
               reactor_state[k],
               reactor_derivative[k],
               exact_derivative[k]); 
      }
      fflush(stdout);
    }
  }
  printf("[%s] PASSED  %d/%d (%7.3f%%)\n",
         argv[0],
         num_passed,
         num_runs,
         100.0*(double)num_passed/(double)num_runs);
  fflush(stdout);

  return 0;
}
static int GetExactDerivative(const double atomic_mass_fraction,
                              const double temperature,
                              const double pressure,
                              std::vector<double> *derivative)
{
  const double gas_constant = zerork::NIST_RU;
  const double mix_molecular_weight = 
    2.0*ATOMIC_MASS/(atomic_mass_fraction+1.0);
  double mix_specific_heat;
  double delta_enthalpy;

  // mixture properties are divided by (Ru/atomic species weight)
  mix_specific_heat = atomic_mass_fraction*ATOMIC_A0 + 
    0.5*(1.0-atomic_mass_fraction)*DIATOMIC_A0;
  delta_enthalpy = (ATOMIC_A0*temperature + ATOMIC_A5) -
    0.5*(DIATOMIC_A0*temperature + DIATOMIC_A5);

  derivative->clear();
  derivative->assign(NUM_STATES,0.0);

  // dy(atomic)/dt
  derivative->at(0) = -4.0*RATE_CONST*pressure/(gas_constant*temperature)*
    atomic_mass_fraction*atomic_mass_fraction/(1.0+atomic_mass_fraction);
  // dy(diatomic)/dt
  derivative->at(1) = -derivative->at(0);

  // d(T/Tref)/dt
  derivative->at(3) = -delta_enthalpy*derivative->at(0)/
    (mix_specific_heat*REF_TEMPERATURE);
  // dv/dt with relative volume v
  derivative->at(2) = 
    gas_constant*REF_TEMPERATURE/(pressure*mix_molecular_weight)*
    derivative->at(3) +
    gas_constant*temperature/pressure*(1.0/ATOMIC_MASS - 0.5/ATOMIC_MASS)*
    derivative->at(0);
    
  return 0.0;
}
static int GetReactorDerivative(const double atomic_mass_fraction,
                                const double temperature,
                                const double pressure,
                                std::vector<double> *state,
                                std::vector<double> *derivative)
{
  double mix_molecular_weight;
  std::vector<double> molecular_weight;

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
  ConstPressureReactor reactor(load_mech.c_str(),
                               load_therm.c_str(),
                               null_filename,
                               COMPRESSED_COL_STORAGE,
                               pressure);
  const double gas_constant = reactor.GetGasConstant();

  reactor.SetReferenceTemperature(REF_TEMPERATURE);
  derivative->clear();
  derivative->assign(NUM_STATES,0.0);

  state->clear();
  state->assign(NUM_STATES, 0.0);

  molecular_weight.clear();
  molecular_weight.assign(NUM_SPECIES, 0.0);
  reactor.GetSpeciesMolecularWeight(&molecular_weight[0]);

  mix_molecular_weight = 1.0/(atomic_mass_fraction/molecular_weight[0] +
                              (1.0-atomic_mass_fraction)/molecular_weight[1]);

  // set the inital state
  state->at(0) = atomic_mass_fraction;
  state->at(1) = 1.0-atomic_mass_fraction;
  state->at(2) = gas_constant*temperature/(pressure*mix_molecular_weight);
  state->at(3) = temperature/REF_TEMPERATURE;

  reactor.GetTimeDerivative(0.0,&state->at(0),&derivative->at(0));
 
  return 0;
}

static int RandomInt(const int a, const int b)
{
  double span = (double)(b-a);
  double random_number = (double)a + zerork::utilities::random01()*span;
  return (int)random_number;
}

static void RandomNormalizedVector(const size_t num_elements,
                                   const size_t num_digits, 
                                   std::vector<double> *v)
{
  double vector_sum = 0.0;
  double rounded_value,remainder;
  double digit_mask = pow(10.0,num_digits);
  v->clear();
  v->assign(num_elements,0.0);
  for(size_t j=0; j<num_elements; ++j) {
    v->at(j) = zerork::utilities::random01();
    vector_sum += v->at(j);
  }
  // normalize
  for(size_t j=0; j<num_elements; ++j) {
    v->at(j) /= vector_sum;
  }
  if(0 < num_digits && num_digits < 9 && num_elements > 0) { 
    // otherwise there will be issues representing the integers
    vector_sum = 0.0;
    for(size_t j=0; j<num_elements-1; ++j) {
      // round down to always add positive numbers
      rounded_value = ((int)(digit_mask*v->at(j)))/digit_mask;
      remainder = v->at(j) - rounded_value;
      v->at(j) = rounded_value;
      v->at(j+1) += remainder;      
      vector_sum += v->at(j);
    }
    // use the last element as the garbage collector to enforce normalization
    v->at(num_elements-1) = 1.0-vector_sum;

    // double check normalization
    vector_sum = 0.0;
    for(size_t j=0; j<num_elements; ++j) {
      vector_sum += v->at(j);
    }
    if(fabs(vector_sum-1.0) > 1.0e-12) {
      RandomNormalizedVector(num_elements,
                             num_digits, 
                             v);
    }   
  }
}

bool SameVector(const std::vector<double>A, 
                const std::vector<double>B,
                const double rtol,
                const double atol)
{
  if(A.size() != B.size()) {
    return false;
  }
  for(size_t j=0; j<A.size(); ++j) {
    double error = fabs(A[j]-B[j]);
    double value = fabs(0.5*(A[j]+B[j]));
    if(error > value*rtol && error > atol) {
      return false;
    }
  }
  return true;
}
