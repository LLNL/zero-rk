#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <string>

#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros
#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.
#if defined SUNDIALS2
#include <cvode/cvode_dense.h>      // prototypes & constants for CVDense
#elif defined SUNDIALS3
#include <cvode/cvode_direct.h>
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#elif defined SUNDIALS4
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#endif

#include <reactor/const_pressure_reactor.h>
#include <utilities/math_utilities.h>
#include <utilities/file_utilities.h>

using zerork::utilities::null_filename;

// TEST CONSTANTS
const double PRESSURE        = 2.5e6;  // [Pa]
const double MIN_TEMPERATURE =  300.0; // [K]
const double MAX_TEMPERATURE = 3500.0; // [K]

const double TEST_ATOL = 1.0e-18;
const double TEST_RTOL = 1.0e-8;

const double CVODE_ATOL = 1.0e-20;
const double CVODE_RTOL = 1.0e-10;
const int MAX_CVODE_STEPS = 100000;
const double MAX_CVODE_TIME = 1.0e3;

const char TEST_MECH[]  = "mechanisms/ideal/hydrogen_recombination.mech";
const char TEST_THERM[] = "mechanisms/ideal/const_specific_heat.therm";

// Thermodynamic property definition taken for hydrogen
const double ATOMIC_A0   =  2.5;
const double ATOMIC_A5   =  2.54737667e+04;
const double DIATOMIC_A0 =  3.5;
const double DIATOMIC_A5 = -1.04352500e+03; 

static double GetFinalTemperature(const double atomic_mass_fraction, 
                                  const double temperature);

static double MarchToFinalTemperature(const double atomic_mass_fraction, 
                                      const double temperature,
                                      const int max_steps,
                                      double *final_time,
                                      int *final_steps);

static int RandomInt(const int a, const int b);

static int ConstPressureRHS(double t, 
                            N_Vector y, 
                            N_Vector ydot,
                            void *user_data);


int main(int argc, char *argv[])
{
  int seed = -1;     // default indicates time(0) used as seed
  int num_runs = 1;
  int num_passed = 0;
  double atomic_mass_fraction;
  double initial_temperature;
  double final_temperature_exact;
  double final_temperature_reactor;
  double error;
  double final_marching_time;
  int    final_marching_steps;
  std::string outcome;

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

  for(int j=0; j<num_runs; ++j) {
    initial_temperature = (double)RandomInt((int)MIN_TEMPERATURE,
                                            (int)MAX_TEMPERATURE);
    atomic_mass_fraction = 1.0e-4*RandomInt(0,10000);
    final_temperature_exact = GetFinalTemperature(atomic_mass_fraction,
                                                  initial_temperature);
    final_temperature_reactor = MarchToFinalTemperature(atomic_mass_fraction,
                                                        initial_temperature,
                                                        MAX_CVODE_STEPS,
                                                        &final_marching_time,
                                                        &final_marching_steps);

    error = final_temperature_reactor - final_temperature_exact;
    if(fabs(error) < TEST_ATOL || 
       fabs(error) < TEST_RTOL*fabs(final_temperature_exact)) {
      outcome = "PASSED";
    } else {
      outcome = "FAILED";
    }


    printf("[%s] RUN %2d: %s T_f(error) %10.3e [K] (final march time %10.3e steps %d)\n",
           argv[0],
           j+1,
           outcome.c_str(),
           error,
           final_marching_time,
           final_marching_steps);
    if(outcome=="FAILED") {
      printf("    RUN %2d: y_0[H] = %6.4f  T_0 = %6.1f [K]\n",
             j+1,
             atomic_mass_fraction,
             initial_temperature);
      printf("    RUN %2d: T_f(exact) = %.18g [K]\n",
             j+1,final_temperature_exact);
      printf("    RUN %2d: T_f(test)  = %.18g [K]\n",
	     j+1,final_temperature_reactor);
    } else {
      ++num_passed;
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

double GetFinalTemperature(const double atomic_mass_fraction, 
                           const double temperature)
{
  //double molecular_weight = 2.0-atomic_mass_fraction; // normalized by the
                                                      // atomic mass
  // enthalpy multiplied by the diatomic molecular weight and 
  // divided by the gas constant, units of [K]
  double enthalpy = 
    2.0*atomic_mass_fraction*(ATOMIC_A0*temperature + ATOMIC_A5) +
    (1.0-atomic_mass_fraction)*(DIATOMIC_A0*temperature + DIATOMIC_A5);

  double final_temperature = (enthalpy-DIATOMIC_A5)/DIATOMIC_A0;
  return final_temperature;
}

static double MarchToFinalTemperature(const double atomic_mass_fraction, 
                                      const double temperature,
                                      const int max_steps,
                                      double *final_time,
                                      int *final_steps)
{
  const double ref_temperature = 1000.0;
  double molecular_weight[2];
  double mix_molecular_weight;
  double current_time = 0;
  double final_temperature;
  int flag = 0;
  int num_steps = 0;
  N_Vector state;
  double *state_data;
  void *cvode_ptr;

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
  ConstPressureReactor user_data(load_mech.c_str(),
                                 load_therm.c_str(),
                                 null_filename,
                                 COMPRESSED_COL_STORAGE,
                                 PRESSURE);
  *final_time  = 0.0;
  *final_steps = 0;

  user_data.SetReferenceTemperature(ref_temperature);
  user_data.GetSpeciesMolecularWeight(&molecular_weight[0]);
  mix_molecular_weight = 1.0/(atomic_mass_fraction/molecular_weight[0]+
                              (1.0-atomic_mass_fraction)/molecular_weight[1]);

  // initialize the state vector
  state = N_VNew_Serial(4);
  state_data = NV_DATA_S(state);
  state_data[0] = atomic_mass_fraction;
  state_data[1] = 1.0-atomic_mass_fraction;
  // relative volume
  state_data[2] = 8314.4598*temperature/(PRESSURE*mix_molecular_weight);
  state_data[3] = temperature/ref_temperature; 
  
  // initialize cvode
#if defined SUNDIALS2 || defined SUNDIALS3
  cvode_ptr = CVodeCreate(CV_BDF, CV_NEWTON);
#elif defined SUNDIALS4
  cvode_ptr = CVodeCreate(CV_BDF);
#endif
  CVodeInit(cvode_ptr, ConstPressureRHS, 0.0, state);
#if defined SUNDIALS2
  CVDense(cvode_ptr, 4);  
#elif defined SUNDIALS3
  SUNMatrix A = SUNDenseMatrix(4, 4);
  SUNLinearSolver LS = SUNDenseLinearSolver(state, A);
  CVDlsSetLinearSolver(cvode_ptr, LS, A);
#elif defined SUNDIALS4
  SUNMatrix A = SUNDenseMatrix(4, 4);
  SUNLinearSolver LS = SUNDenseLinearSolver(state, A);
  CVodeSetLinearSolver(cvode_ptr, LS, A);
#endif
  CVodeSStolerances(cvode_ptr, CVODE_RTOL, CVODE_ATOL);
  CVodeSetUserData(cvode_ptr, (void *)&user_data);
  CVodeSetInitStep(cvode_ptr, CVODE_ATOL); 
 
  while(num_steps < max_steps && state_data[0] > CVODE_ATOL && flag >= 0) {
    flag = CVode(cvode_ptr, 
                 MAX_CVODE_TIME, 
                 state, 
                 &current_time,
		 CV_ONE_STEP);
    ++num_steps;
    //mix_molecular_weight = 1.0/(state_data[0]/molecular_weight[0] +
    //                            state_data[1]/molecular_weight[1]);
    //double current_pressure = 
    //  8314.4598*state_data[3]*ref_temperature/
    //  (state_data[2]*mix_molecular_weight);
    //printf("step %4d: %14.7e  %14.7e  %14.7e  %14.7e  %14.7e %14.7e\n",
    //       num_steps,
    //       state_data[0],
    //       state_data[1],
    //       state_data[2],
    //       state_data[3]*ref_temperature,
    //       state_data[0]+state_data[1]-1.0, // mass fraction error
    //       current_pressure-PRESSURE);      // const pressure error
    //fflush(stdout);
            
  }
  if(flag < 0 || num_steps >= max_steps) {
    printf("WARNING: In MarchToFinalTemperature(...),\n");
    printf("         integration stopped at\n");
    printf("                         time = %.18g [s]\n",current_time);
    printf("                        steps = %d\n", num_steps);
    printf("         atomic mass fraction = %.18g\n",state_data[0]);
    printf("                  temperature = %.18g [K]\n",
           state_data[3]*ref_temperature);
    fflush(stdout);
  }

  final_temperature = state_data[3]*ref_temperature;
  *final_time  = current_time;  
  *final_steps = num_steps;

  N_VDestroy_Serial(state);
  CVodeFree(&cvode_ptr);  
#if defined SUNDIALS3 || defined SUNDIALS4
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
#endif
  return final_temperature;
}



static int RandomInt(const int a, const int b)
{
  double span = (double)(b-a);
  double random_number = (double)a + zerork::utilities::random01()*span;
  return (int)random_number;
}

static int ConstPressureRHS(double t, 
                            N_Vector y, 
                            N_Vector ydot,
                            void *user_data)
{
  ConstPressureReactor *reactor = (ConstPressureReactor *)user_data;
  reactor->GetTimeDerivative(t,NV_DATA_S(y),NV_DATA_S(ydot));
  return 0;
}

