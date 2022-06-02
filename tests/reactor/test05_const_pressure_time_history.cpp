#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cstdlib>

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

#include <zerork/constants.h>
#include <reactor/const_pressure_reactor.h>
#include <utilities/math_utilities.h>
#include <utilities/file_utilities.h>

using zerork::utilities::null_filename;

typedef struct
{
  double atomic_mass_fraction_;
  double temperature_;
  double pressure_;
} InitialConditions;

typedef struct
{
  double rate_constant_;
  double atomic_a0_;
  double atomic_a5_;
  double diatomic_a0_;
  double diatomic_a5_;
} MechanismParameters;


// TEST CONSTANTS
const double MIN_PRESSURE    = 1.0e5;  // [Pa]
const double MAX_PRESSURE    = 1.0e7;
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
const double RATE_CONSTANT = 2.0e9;
const double ATOMIC_A0   =  2.5;
const double ATOMIC_A5   =  2.54737667e+04;
const double DIATOMIC_A0 =  3.5;
const double DIATOMIC_A5 = -1.04352500e+03; 


static double GetMassFractionIntegral(const InitialConditions &ic,
                                      const MechanismParameters &mech,
                                      const double atomic_mass_fraction);

static double MarchToTime(const InitialConditions &ic,
                          const MechanismParameters &mech,
                          const double time_target,
                          const int max_steps);

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
  double initial_pressure;
  double reduction;
  double error;
  InitialConditions ic;
  MechanismParameters mech;
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

  mech.rate_constant_ = RATE_CONSTANT;
  mech.atomic_a0_     = ATOMIC_A0;
  mech.atomic_a5_     = ATOMIC_A5;
  mech.diatomic_a0_   = DIATOMIC_A0;
  mech.diatomic_a5_   = DIATOMIC_A5;
 
  for(int j=0; j<num_runs; ++j) {
    initial_temperature = (double)RandomInt((int)MIN_TEMPERATURE,
                                            (int)MAX_TEMPERATURE);
    initial_pressure = (double)RandomInt((int)MIN_PRESSURE,
                                         (int)MAX_PRESSURE);
    atomic_mass_fraction = 1.0e-4*RandomInt(0,10000);
    reduction = pow(10.0,-RandomInt(1,8));

    ic.temperature_          = initial_temperature;
    ic.pressure_             = initial_pressure;
    ic.atomic_mass_fraction_ = atomic_mass_fraction;

    double F0 = GetMassFractionIntegral(ic, mech, atomic_mass_fraction);
    double F1 = GetMassFractionIntegral(ic, mech, reduction*atomic_mass_fraction);
    double elapsed_time = F1-F0;

    double final_atomic_mass_fraction = MarchToTime(ic,
                                                    mech,
                                                    elapsed_time,
                                                    MAX_CVODE_STEPS);

    error = final_atomic_mass_fraction - reduction*atomic_mass_fraction;

    if(fabs(error) < fabs(TEST_ATOL) || 
       fabs(error) < fabs(TEST_RTOL*reduction*atomic_mass_fraction)) {

      outcome = "PASSED";
    } else {
      outcome = "FAILED";
    }


    printf("[%s] RUN %2d: %s  t=%8.2e [s]  y_f(error)=%9.2e  y_f(rel. error)=%9.2e\n",
           argv[0],
           j+1,
           outcome.c_str(),
           elapsed_time,
           error,
           fabs(error/(reduction*atomic_mass_fraction)));
    if(outcome=="FAILED") {
      printf("    RUN %2d: y_0[H] = %6.4f  T_0 = %6.1f [K]  p_0 = %4.1f [Pa]\n",
             j+1,
             atomic_mass_fraction,
             initial_temperature,
             initial_pressure);
      printf("    RUN %2d: y_f(exact) = %.18g [-]\n",
             j+1,reduction*atomic_mass_fraction);
      printf("    RUN %2d: y_f(test)  = %.18g [-]\n",
    	     j+1,final_atomic_mass_fraction);
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

// double GetFinalTemperature(const double atomic_mass_fraction, 
//                            const double temperature)
// {
//   //double molecular_weight = 2.0-atomic_mass_fraction; // normalized by the
//                                                       // atomic mass
//   // enthalpy multiplied by the diatomic molecular weight and 
//   // divided by the gas constant, units of [K]
//   double enthalpy = 
//     2.0*atomic_mass_fraction*(ATOMIC_A0*temperature + ATOMIC_A5) +
//     (1.0-atomic_mass_fraction)*(DIATOMIC_A0*temperature + DIATOMIC_A5);

//   double final_temperature = (enthalpy-DIATOMIC_A5)/DIATOMIC_A0;
//   return final_temperature;
// }

static double MarchToTime(const InitialConditions &ic,
                          const MechanismParameters &mech,
                          const double time_target,
                          const int max_steps)
{
  const double ref_temperature = 1000.0;
  double molecular_weight[2];
  double mix_molecular_weight;
  double current_time;
  int flag = 0;
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
                                 ic.pressure_);

  user_data.SetReferenceTemperature(ref_temperature);
  user_data.GetSpeciesMolecularWeight(&molecular_weight[0]);
  mix_molecular_weight = 
    1.0/(ic.atomic_mass_fraction_/molecular_weight[0]+
         (1.0-ic.atomic_mass_fraction_)/molecular_weight[1]);

  // initialize the state vector
  current_time = 0.0;
  state = N_VNew_Serial(4);
  state_data = NV_DATA_S(state);
  state_data[0] = ic.atomic_mass_fraction_;
  state_data[1] = 1.0-ic.atomic_mass_fraction_;
  // relative volume
  state_data[2] = 
    user_data.GetGasConstant()*ic.temperature_/
    (ic.pressure_*mix_molecular_weight);
  state_data[3] = ic.temperature_/ref_temperature; 
  
  // initialize cvode
#if defined SUNDIALS2 || defined SUNDIALS3
  cvode_ptr = CVodeCreate(CV_BDF, CV_NEWTON);
#elif defined SUNDIALS4
  cvode_ptr = CVodeCreate(CV_BDF);
#endif
  CVodeInit(cvode_ptr, ConstPressureRHS, current_time, state);
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
  CVodeSetMaxNumSteps(cvode_ptr, max_steps); 

  flag = CVode(cvode_ptr, 
               time_target, 
               state, 
               &current_time, 
               CV_NORMAL);
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
            

  if(flag < 0) {
    printf("WARNING: In MarchToFinalTemperature(...),\n");
    printf("         integration stopped at\n");
    printf("                         time = %.18g [s]\n",current_time);
    printf("         atomic mass fraction = %.18g\n",state_data[0]);
    printf("                  temperature = %.18g [K]\n",
           state_data[3]*ref_temperature);
    fflush(stdout);
  }

  double final_result = state_data[0];

  N_VDestroy_Serial(state);
  CVodeFree(&cvode_ptr);  
#if defined SUNDIALS3 || defined SUNDIALS4
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
#endif
  return final_result;
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

// Solves the mass fraction integral
//
// F(y) = Integral[ C*(m1*y + b1)/(m2*y + b2) * (y+1)/y^2 ]dy
// 
// where y is the atomic mass fraction
//
// F(y(t)) - F(y(0)) = t
static double GetMassFractionIntegral(const InitialConditions &ic,
                                      const MechanismParameters &mech,
                                      const double atomic_mass_fraction)
{
  const double C = -zerork::NIST_RU/(4.0*mech.rate_constant_*ic.pressure_);
  const double m1 = -mech.atomic_a5_ + 0.5*mech.diatomic_a5_;
  const double b1 = ic.atomic_mass_fraction_*
    (mech.atomic_a0_*ic.temperature_ + mech.atomic_a5_) +
    0.5*(1.0-ic.atomic_mass_fraction_)*
    (mech.diatomic_a0_*ic.temperature_ + mech.diatomic_a5_) - 
    0.5*mech.diatomic_a5_;
  const double m2 = mech.atomic_a0_ - 0.5*mech.diatomic_a0_;
  const double b2 = 0.5*mech.diatomic_a0_;
  const double y = atomic_mass_fraction; // simplifying formulas

  double result = C/(b2*b2*m2*y)*
    (-b1*b2*m2 + (b2*m1 + b1*(b2-m2))*m2*y*log(y) +
     (b2-m2)*(b2*m1 - b1*m2)*y*log(m2*y + b2));
 
  return result;
}
