#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <nvector/nvector_serial.h> // serial N_Vector types, fcts., and macros
#include <cvode/cvode.h>            // prototypes for CVODE fcts. and consts.

#include <ext/cvode_user_sparselu.h>

const int NUM_STATES = 3;
const int NUM_NON_ZEROS = 5;
const int test_row_index[]={0,0,1,1,2};
const int test_col_sum[]={0,1,3,5};
const double test_rtol = 1.0e-12;
const double test_atol = 1.0e-20;
const double test_tmax = 1.0;
const double test_dt   = 0.02;

int sparse_ode_rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int sparse_ode_jac(int N,
                   realtype t,
                   N_Vector y,
                   N_Vector fy,
                   double user_jac[],
                   void *user_data,
                   N_Vector tmp1,
                   N_Vector tmp2,
                   N_Vector tmp3);

void exact_solution(const double t, const double y_init[], double y[]);


int main(int argc, char *argv[])
{
  int j;
  int flag;
  int num_jac,num_factor,num_solve;
  double y_init[NUM_STATES];
  double y_exact[NUM_STATES];
  double tnext = test_dt;
  double tcurr = 0.0;
  N_Vector system_state;
  void *cvode_mem;

  // Set the vector of initial values
  system_state = N_VNew_Serial(NUM_STATES);
  for(j=0; j<NUM_STATES; ++j) {
    y_init[j] = NV_Ith_S(system_state,j) = (double)(j+1);
  }

  // Create a CVode memory object 
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  // Initialize the CVode solver
  flag = CVodeInit(cvode_mem, sparse_ode_rhs, 0.0, system_state);

  // Specify the integration tolerances
  flag = CVodeSStolerances(cvode_mem, test_rtol, test_atol);

  
  // Specify optional inputs using CVodeSet*
  // :
  // :

  // Attach User defined SuperLU linear solver
  flag = CVUserSuperLU(cvode_mem,
                       NUM_STATES,
                       NUM_NON_ZEROS,
                       test_row_index,
                       test_col_sum,
                       sparse_ode_jac);
  if(flag != CV_SUCCESS) {
    printf("ERROR: CVUserSuperLU(...) returned flag = %d\n",flag);
    exit(-1);
  }
  
  printf("# time [s]          y[0]                y[1]                y[2]  rel err y[0]  rel err y[1]  rel err y[2]\n");
  while(tcurr < test_tmax - 0.5*test_dt) {
    // March solution
    flag = CVode(cvode_mem,
                 tnext, 
                 system_state,
                 &tcurr,
                 CV_NORMAL);
    if(flag != CV_SUCCESS) {
      printf("ERROR: CVode(...) returned flag = %d\n",flag);
      exit(-1);
    }
    exact_solution(tcurr, y_init, y_exact);
    printf("%4.2f  %18.12e  %18.12e  %18.12e  %12.5e  %12.5e  %12.5e\n",
           tcurr,   
           NV_Ith_S(system_state,0),
           NV_Ith_S(system_state,1),
           NV_Ith_S(system_state,2),
           (NV_Ith_S(system_state,0)-y_exact[0])/y_exact[0],
           (NV_Ith_S(system_state,1)-y_exact[1])/y_exact[1],
	   (NV_Ith_S(system_state,2)-y_exact[2])/y_exact[2]);
    fflush(stdout);
    tnext+=test_dt;
  }
  CVUserSuperLUGetNumJacEvals(  cvode_mem, &num_jac);
  CVUserSuperLUGetNumJacFactors(cvode_mem, &num_factor);
  CVUserSuperLUGetNumJacSolves( cvode_mem, &num_solve);
  

  printf("# Summary:\n");
  printf("# Number of CVUserSuperLU jacobian evaluations    : %d\n",
         num_jac);
  printf("# Number of CVUserSuperLU jacobian factorizations : %d\n",
         num_factor);
  printf("# Number of CVUserSuperLU jacobian back solves    : %d\n",
         num_solve);
  

  N_VDestroy_Serial(system_state);
  CVodeFree(&cvode_mem);

  return 0;
}

int sparse_ode_rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  double y0 = NV_Ith_S(y,0);
  double y1 = NV_Ith_S(y,1);
  double y2 = NV_Ith_S(y,2);

  NV_Ith_S(ydot,0) =  4.0*y0 - y1;
  NV_Ith_S(ydot,1) = -2.0*y1 + 3.0*y2;
  NV_Ith_S(ydot,2) =  -5.0*y2;
  return 0;
}
int sparse_ode_jac(int N,
                   realtype t,
                   N_Vector y,
                   N_Vector fy,
                   double user_jac[],
                   void *user_data,
                   N_Vector tmp1,
                   N_Vector tmp2,
                   N_Vector tmp3)
{
  user_jac[0] =  4.0; // J[0][0] =  4.0
  user_jac[1] = -1.0; // J[0][1] = -1.0
  user_jac[2] = -2.0; // J[1][1] = -2.0
  user_jac[3] =  3.0; // J[1][2] =  3.0
  user_jac[4] = -5.0; // J[2][2] = -5.0
  return 0;
}

void exact_solution(const double t, const double y_init[], double y[])
{
  double z_init[NUM_STATES];
  double z[NUM_STATES];

  z_init[0] = y_init[0] - y_init[1]/6.0 - y_init[2]/18.0;
  z_init[1] =             y_init[1]     + y_init[2];
  z_init[2] =                             y_init[2];

  z[0] = z_init[0]*exp( 4.0*t);
  z[1] = z_init[1]*exp(-2.0*t);
  z[2] = z_init[2]*exp(-5.0*t);
 
  y[0] = z[0] + z[1]/6.0 - z[2]/9.0;
  y[1] =        z[1]     - z[2];
  y[2] =                   z[2];
}
