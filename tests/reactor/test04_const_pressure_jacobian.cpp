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

#include "test04_utilities.h"

// TEST CONSTANTS
const double MIN_PRESSURE    = 1.0e5;  // [Pa]
const double MAX_PRESSURE    = 1.0e7;  // [Pa]

const double MIN_TEMPERATURE =  300.0; // [K]
const double MAX_TEMPERATURE = 3500.0; // [K]

const double REF_TEMPERATURE = 3000.0; // [K]

const double TEST_RTOL = 1.0e-6;
const double TEST_ATOL = 1.0e-12;

const double PERTURB_RDELTA = 1.0e-8;
const double PERTURB_ADELTA = 1.0e-14;

const char TEST_MECH[]  = "mechanisms/ideal/hydrogen_no_falloff.mech";
const char TEST_THERM[] = "mechanisms/ideal/const_specific_heat.therm";
//const char TEST_MECH[]  = "mechanisms/hydrogen/h2_v1b_mech.txt";
//const char TEST_THERM[] = "mechanisms/hydrogen/h2_v1a_therm.txt";
//const char TEST_MECH[]  = "mechanisms/ideal/hydrogen_recombination.mech";
//const char TEST_THERM[] = "mechanisms/ideal/const_specific_heat.therm";


static int GetDirectJacobian(ConstPressureReactor &reactor,
                             const std::vector<double> &state,
                             std::vector<double> *matrix);

static int RandomInt(const int a, const int b);

static void RandomNormalizedVector(const size_t num_elements,
                                   const size_t num_digits, 
                                   std::vector<double> *v);

double GetLargestWeightedError(const std::vector<double>A, 
                               const std::vector<double>B,
                               const double rtol,
                               const double atol,
                               int *largest_id,
                               std::vector<double> *weighted_error);

double GetWeightedError(const std::vector<double>A, 
                        const std::vector<double>B,
                        const double rtol,
                        const double atol,
                        const std::vector<double>roundoff,
                        std::vector<double> *weighted_error);

double GetL2Norm(const std::vector<double>A);
double GetL2Difference(const std::vector<double>A, 
                       const std::vector<double>B);
double GetL2DifferenceIgnoreEpsilon(const std::vector<double>A, 
                                    const std::vector<double>B,
                                    const std::vector<double>epsilon);

bool SameVector(const std::vector<double>A, 
                const std::vector<double>B,
                const double rtol,
                const double atol);

bool SameVector(const std::vector<double>A, 
                const std::vector<double>B,
                const double rtol,
                const double atol,
                const std::vector<double>roundoff);



int ConstPressureDerivative(const size_t num_equations,
                            const double t, 
                            const double y[],
                            double f[],
                            void *params);


int main(int argc, char *argv[])
{
  int seed = -1;     // default indicates time(0) used as seed
  int num_runs = 1;
  int num_passed = 0;
  std::string outcome;
  //double initial_temperature, initial_pressure;
  std::vector<double> mass_fraction;
  std::vector<double> initial_state, initial_derivative;
  std::vector<double> molecular_weight;
  std::vector<double> matrix_from_direct, matrix_from_rhs;
  std::vector<double> roundoff_error, truncation_error;
  double mass_sum, mix_molecular_weight;
  ConstPressureReactor *user_data;

  if(argc >= 3) {
    seed = atoi(argv[2]);
  } 
  if(argc >= 2) {
    num_runs = atoi(argv[1]);
  }
  if(num_runs <= 0) {
    num_runs = 1;
  }
  if(seed == -1) {
    seed = time(0);
  }
  zerork::utilities::random01seed(seed);

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

  for(int j=0; j<num_runs; ++j) {
    double initial_temperature = (double)RandomInt((int)MIN_TEMPERATURE,
                                                   (int)MAX_TEMPERATURE);
    double initial_pressure = (double)RandomInt((int)MIN_PRESSURE,
                                                (int)MAX_PRESSURE);
    user_data = new ConstPressureReactor(load_mech.c_str(),
                                         load_therm.c_str(),
                                         null_filename,
                                         COMPRESSED_COL_STORAGE,
                                         initial_pressure);
    user_data->SetReferenceTemperature(REF_TEMPERATURE);

    const int num_species = user_data->GetNumSpecies();
    const int num_states  = user_data->GetNumStates();
    RandomNormalizedVector(num_species, 
                           3, // number of decimal places
                           &mass_fraction);

    initial_state.clear();
    initial_state.assign(num_states, 0.0);
    initial_derivative.clear();
    initial_derivative.assign(num_states, 0.0);
 
    molecular_weight.clear();
    molecular_weight.assign(num_species, 0.0);
    user_data->GetSpeciesMolecularWeight(&molecular_weight[0]);

    mass_sum = 0.0;
    for(int k=0; k<num_species; ++k) {
      initial_state[k] = mass_fraction[k];
      mass_sum += mass_fraction[k]/molecular_weight[k];
    }
    mix_molecular_weight = 1.0/mass_sum;        
    initial_state[num_species] = user_data->GetGasConstant()*
      initial_temperature/(mix_molecular_weight*initial_pressure);
    initial_state[num_species+1] = 
      initial_temperature/user_data->GetReferenceTemperature();

    GetDirectJacobian(*user_data,
                      initial_state,
                      &matrix_from_direct);
    //GetJacobianFromRHS(*user_data,
    //                   initial_state,
    //                   &matrix_from_rhs,
    //                   &roundoff_error);

    matrix_from_rhs.clear();
    matrix_from_rhs.assign(num_states*num_states,0.0);

    truncation_error.clear();
    truncation_error.assign(num_states*num_states,0.0);

    roundoff_error.clear();
    roundoff_error.assign(num_states*num_states,0.0);

    int error_flag = GetJacobian(//FORWARD_FIRST_ORDER,
                                 CENTRAL_FOURTH_ORDER,
                                 ConstPressureDerivative,
                                 PERTURB_RDELTA,
                                 PERTURB_ADELTA,
                                 (size_t)num_states,
                                 0.0,
                                 &initial_state[0],
                                 user_data,
                                 &matrix_from_rhs[0],
                                 &truncation_error[0],
                                 &roundoff_error[0]);
    if(error_flag != 0) {
      printf("WARNING: GetJacobian(...) returned error flag = %d\n",
             error_flag);
      fflush(stdout);
    }


    double frobenius_norm  = GetL2Norm(matrix_from_direct);
    double difference_norm = GetL2Difference(matrix_from_direct,
                                             matrix_from_rhs);
    double difference_norm_ignore_roundoff = 
      GetL2DifferenceIgnoreEpsilon(matrix_from_direct,
                                   matrix_from_rhs,
                                   roundoff_error);


    if(difference_norm_ignore_roundoff <= TEST_ATOL ||
       difference_norm_ignore_roundoff <= frobenius_norm*TEST_RTOL ) {
      outcome = "PASSED";
      ++num_passed;
    } else {
      outcome = "FAILED";
    }
    printf("[%s] RUN %2d: %s   Frobenius norms |J| = %10.3e |error| = %10.3e |error| = %10.3e (exclude roundoff)",
           argv[0],
           j+1,
           outcome.c_str(),
           frobenius_norm,
           difference_norm,
           difference_norm_ignore_roundoff);
    if(outcome == "FAILED") {
      std::vector<double> weighted_error;
      int max_id;
      double max_error = GetLargestWeightedError(matrix_from_direct,
                                                 matrix_from_rhs,
                                                 TEST_RTOL,
                                                 TEST_ATOL,
                                                 &max_id,
                                                 &weighted_error);
      printf("  max error = %9.2e at (%d, %d)",
             max_error, 1+max_id%num_states,1+max_id/num_states);

      max_error = GetWeightedError(matrix_from_direct,
                                   matrix_from_rhs,
                                   TEST_RTOL,
                                   TEST_ATOL,
                                   roundoff_error,
                                   &weighted_error);
      printf(" exclude roundoff %9.2e\n", max_error);

      user_data->GetTimeDerivative(0.0,
                                   &initial_state[0],
                                   &initial_derivative[0]);
                                   
      printf("  Initial State:\n");
      for(int k=0; k<num_states; ++k) {
        printf("    y[%2d]: %24.16e dy/dt = %24.16e (%s)\n",
               k+1,
               initial_state[k],
               initial_derivative[k],
               user_data->GetNameOfStateId(k)); 
      }

      printf("  Jacobian Comparison:\n");
      for(int k=0; k<num_states; ++k) {   // column-k
        for(int m=0; m<num_states; ++m) { // row-m
          printf("    J[%2d,%2d]: %10.3e %10.3e %10.3e  %24.16e  %24.16e  (%s,%s)\n",
                 m+1,
                 k+1,
                 weighted_error[m+num_states*k],
                 roundoff_error[m+num_states*k],
                 fabs(matrix_from_direct[m+num_states*k]-
                      matrix_from_rhs[m+num_states*k]),
                 matrix_from_direct[m+num_states*k],
                 matrix_from_rhs[m+num_states*k],
                 user_data->GetNameOfStateId(m),
                 user_data->GetNameOfStateId(k));
        }
      }
    }  else {
      printf("\n");
    } 
    delete user_data;
  }
  printf("[%s] PASSED  %d/%d (%7.3f%%)\n",
         argv[0],
         num_passed,
         num_runs,
         100.0*(double)num_passed/(double)num_runs);
  fflush(stdout);

  return 0;
}

static int GetDirectJacobian(ConstPressureReactor &reactor,
                             const std::vector<double> &state,
                             std::vector<double> *matrix)
{
  const int num_states = reactor.GetNumStates();
  const int num_nonzeros = reactor.GetJacobianSize();
  std::vector<double> sparse_jacobian;
  std::vector<int> row_id, column_id;

  sparse_jacobian.clear();
  sparse_jacobian.assign(num_nonzeros,0.0);
  row_id.clear();
  row_id.assign(num_nonzeros,0);
  column_id.clear();
  column_id.assign(num_nonzeros,0);

  reactor.GetJacobianPattern(&row_id[0],&column_id[0]);
  reactor.GetJacobian(0.0,
                      &state[0],
                      &sparse_jacobian[0]);

  matrix->clear();
  matrix->assign(num_states*num_states,0.0);

  for(int j=0; j<num_nonzeros; ++j) {
    matrix->at(row_id[j]+column_id[j]*num_states) = sparse_jacobian[j];
  }

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

bool SameVector(const std::vector<double>A, 
                const std::vector<double>B,
                const double rtol,
                const double atol,
                const std::vector<double>roundoff)
{
  if(A.size() != B.size()) {
    return false;
  }
  for(size_t j=0; j<A.size(); ++j) {
    double error = fabs(A[j]-B[j]);
    double value = fabs(0.5*(A[j]+B[j]));

    if(error > value*rtol && error > atol && error > fabs(roundoff[j])) {
      // the error must fail both the rtol and atol checks and be larger
      // than the roundoff error
      
      return false;
    }
  }
  return true;

}


//
double GetLargestWeightedError(const std::vector<double>A, 
                               const std::vector<double>B,
                               const double rtol,
                               const double atol,
                               int *largest_id,
                               std::vector<double> *weighted_error)
{
  double current_error;
  double max_error = -1.0e300;
  *largest_id = 0;  
  weighted_error->clear();

  if(A.size() != B.size()) {
    *largest_id = -1;
    return max_error;
  }
  weighted_error->assign(A.size(),0.0);

  for(size_t j=0; j<A.size(); ++j) {
    double value = fabs(0.5*(A[j]+B[j]));
    // cvode error weight
    current_error = fabs(A[j]-B[j])/fabs(rtol*value+atol);
    weighted_error->at(j) = current_error;    

    if(current_error > max_error) {
      max_error = current_error;
      *largest_id = j;
    }
  }
  return max_error;
}

double GetWeightedError(const std::vector<double>A, 
                        const std::vector<double>B,
                        const double rtol,
                        const double atol,
                        const std::vector<double>roundoff,
                        std::vector<double> *weighted_error)
{
  const int min_size = ((A.size() < B.size()) ? A.size() : B.size());
  double max_error=-1.0e300;
  double current_error, current_magnitude;

  weighted_error->clear();
  weighted_error->assign(min_size,0.0);  
  
  for(int j=0; j<min_size; ++j) {
    current_error     = fabs(A[j]-B[j]);
    current_magnitude = 0.5*(A[j]+B[j]); 
    if(current_error > fabs(roundoff[j])) {
      weighted_error->at(j) = 
        current_error/(fabs(rtol*current_magnitude)+fabs(atol));
    }
    if(weighted_error->at(j) > max_error) {
      max_error = weighted_error->at(j);
    }
  }  
  return max_error;
}

double GetL2Norm(const std::vector<double>A)
{
  double sum=0.0;
  for(size_t j=0; j<A.size(); ++j) {
    sum+=A[j]*A[j];
  } 
  return sqrt(sum);
}

double GetL2Difference(const std::vector<double>A, 
                       const std::vector<double>B)
{
  if(A.size() != B.size()) {
    printf("WARNING: In GetL2Difference(A,B),\n");
    printf("         vector lengths are not equal.\n");
    printf("             A.size() = %lu\n",A.size());
    printf("             B.size() = %lu\n",B.size());
    printf("         Returning zero.\n");
    fflush(stdout);
    return 0.0;
  }
  std::vector<double> difference;
  difference.assign(A.size(),0.0);
  for(size_t j=0; j<A.size(); ++j) {
    difference[j] = A[j]-B[j];
  }
  return GetL2Norm(difference);
}
double GetL2DifferenceIgnoreEpsilon(const std::vector<double>A, 
                                    const std::vector<double>B,
                                    const std::vector<double>epsilon)
{
  if(A.size() != B.size()) {
    printf("WARNING: In GetL2Difference(A,B),\n");
    printf("         vector lengths are not equal.\n");
    printf("             A.size() = %lu\n",A.size());
    printf("             B.size() = %lu\n",B.size());
    printf("         Returning zero.\n");
    fflush(stdout);
    return 0.0;
  }
  std::vector<double> difference;
  difference.assign(A.size(),0.0);
  for(size_t j=0; j<A.size(); ++j) {
    difference[j] = A[j]-B[j];
    if(fabs(difference[j]) < fabs(epsilon[j])) {
      difference[j] = 0.0;
    }
  }
  return GetL2Norm(difference);
}


int ConstPressureDerivative(const size_t num_equations,
                            const double t, 
                            const double y[],
                            double f[],
                            void *params)
{
  ConstPressureReactor *reactor = (ConstPressureReactor *)params;
  if((int)num_equations != reactor->GetNumStates()) {
    return -1;
  }
  reactor->GetTimeDerivative(t,y,f);
  return 0;
}
