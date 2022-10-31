#ifndef SEULEX_H
#define SEULEX_H

#include <vector>
#include "sundials/sundials_nvector.h"

namespace seulex_cpp {

typedef int (*seul_deriv_fcn)(double x, N_Vector y, N_Vector dy, void* user_data);
typedef int (*seul_jac_fcn)(double x, N_Vector y, N_Vector fy, void* user_data);
typedef int (*seul_jac_decomp_fcn)(double hji, void* user_data);
typedef int (*seul_jac_solve_fcn)(double x, N_Vector y, N_Vector fy,
                                  N_Vector r, N_Vector z, void* user_data);
typedef int (*seul_output_fcn)(int nsteps, double x, double h, N_Vector y, N_Vector dy, void* user_data);

class seulex {
 public:
  seulex(int, N_Vector);
  virtual ~seulex();

  void set_nmax(int);
  void set_km(int);
  void set_nsequ(int);
  void set_lambda(int);
  void set_nrdens(int);
  void set_uround(double);
  void set_hinit(double);
  void set_hmax(double);
  void set_thet(double);
  void set_step_size_params(double,double);
  void set_order_params(double,double);
  void set_safe1safe2(double,double);
  void set_work_params(double,double,double,double);
  void set_autonomous(bool);
  void set_implicit(bool);
  void set_bandwidth(bool);
  void set_tolerances(N_Vector, N_Vector);
  void set_tolerances(double, double);
  void set_deriv_fcn(seul_deriv_fcn);
  void set_jac_fcn(seul_jac_fcn);
  void set_jac_decomp_fcn(seul_jac_decomp_fcn);
  void set_jac_solve_fcn(seul_jac_solve_fcn);
  void set_output_fcn(seul_output_fcn, void*);
  void set_user_data(void*);
  int solve(double* x_in, double xend, N_Vector y);
  void get_integrator_stats(int* nfcn,int* njac,int* nstep,
                            int* naccpt,int* nrejct,int* ndec,int* nsol);
 private:
  int n;
  int nfcn;
  int njac;
  int nstep;
  int naccept;
  int nreject;
  int ndec;
  int nsol;
// -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
  int nmax;
// -------- KM     MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION
  int km;
// -------- NSEQU     CHOICE OF STEP SIZE SEQUENCE
  int nsequ;
// -------- LAMBDA   PARAMETER FOR DENSE OUTPUT
  int lambda;
// -------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS
  int nrdens;
// -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0
  double uround;
// -------- MAXIMAL STEP SIZE
  double h;
  double hmax;
// ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
  double thet;
  double theta;

// -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION
  double fac1;
  double fac2;
// -------  FAC3, FAC4   PARAMETERS FOR THE ORDER SELECTION
  double fac3;
  double fac4;
// ------- SAFE1, SAFE2 SAFETY FACTORS FOR STEP SIZE PREDICTION
  double safe1;
  double safe2;
// ------- WKFCN,WKJAC,WKDEC,WKSOL  ESTIMATED WORK FOR  FCN,JAC,DEC,SOL
  double wkfcn;
  double wkjac;
  double wkdec;
  double wksol;
  double wkrow;

  N_Vector rtol;
  N_Vector atol;

  bool caljac;
  bool last;
  bool reject;
  bool atov;
  double err;
  double errold;
  double x;

  void check_tolerances(void);
  void create_work_vectors(N_Vector y);
  void destroy_work_vectors();
  void setup_step_sequence(std::vector<int>& nj, std::vector<double>& a);
  int seul(double x, N_Vector y, int jj,
           std::vector<double>& hh, std::vector<double>& w,
           std::vector<int>& nj, std::vector<double>& a);

// Work arrays
  N_Vector yh;
  N_Vector dy;
  N_Vector dyh;
  N_Vector del;
  N_Vector wh;
  N_Vector scal;
  N_Vector tmp1;
  std::vector<N_Vector> t;

  seul_deriv_fcn deriv_fcn;
  seul_jac_fcn jac_fcn;
  seul_jac_decomp_fcn jac_decomp_fcn;
  seul_jac_solve_fcn jac_solve_fcn;
  seul_output_fcn output_fcn;

  void* output_fcn_data;
  void* user_data;
};

} //end namespace seulex_cpp

#endif
