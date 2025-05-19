#ifndef REACTOR_BASE_H_
#define REACTOR_BASE_H_

#include <vector>

#include "optionable.h"

#include "sundials/sundials_nvector.h"
#include "sundials/sundials_direct.h" 
#if defined SUNDIALS3 || defined SUNDIALS4
#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>
#endif

class ReactorBase : public Optionable {
 public:
  ReactorBase() {};
  virtual ~ReactorBase() {};

  virtual void InitializeState(const double reactor_time,
                               const int n_reactors,
                               const double *T,
                               const double *P,
                               const double *mf,
                               const double *dpdt,
                               const double *e_src,
                               const double *y_src) = 0;

  virtual N_Vector& GetStateNVectorRef() = 0;
  virtual void SetBatchMaskNVector(int reactor_idx, N_Vector batch_mask) = 0;

  virtual void GetState(const double reactor_time,
                        double *T,
                        double *P,
                        double *mf) = 0;

  virtual std::vector<double>& GetReactorWeightsRef() = 0;

  virtual void GetAbsoluteToleranceCorrection(N_Vector correction) = 0;

  virtual int GetTimeDerivative(const double reactor_time,
                                N_Vector state,
                                N_Vector derivative) = 0;

  int GetID();
  void SetID(int id);
  virtual void SetStepLimiter(double value) = 0;

#ifdef SUNDIALS2
  virtual int GetJacobianDense(long int N, double t, N_Vector y, N_Vector fy,
                               DlsMat Jac) = 0;
#elif defined SUNDIALS3 || defined SUNDIALS4
  virtual int GetJacobianDense(double t, N_Vector y, N_Vector fy,
                               SUNMatrix Jac) = 0;
#else
#error "Unsupported SUNDIALS version"
#endif

  virtual int GetJacobianDenseRaw(long int N, double t, N_Vector y, N_Vector fy,
                                  double* Jac) = 0;

  virtual int JacobianSetup(double t, N_Vector y, N_Vector fy) = 0;

  virtual int JacobianFactor(double gamma) = 0;

  virtual int JacobianSolve(double t, N_Vector y, N_Vector fy,
                            N_Vector r, N_Vector z) = 0;

  virtual int ComplexJacobianFactor(int k, double alpha, double beta) = 0;

  virtual int ComplexJacobianSolve(int k, N_Vector ax, N_Vector bx) = 0;

  virtual int RootFunction(double t, N_Vector y, double *root_function) = 0;

  virtual int SetRootTime(double t) = 0;
  virtual double GetRootTime() = 0;

  virtual int GetNumStateVariables() = 0;

  virtual int GetNumRootFunctions() = 0;

  virtual int GetNumBatchReactors() = 0;
  virtual int GetMinBatchReactors() = 0;
  virtual int GetMaxBatchReactors() = 0;

  virtual void SetSolveTemperature(bool value) = 0;

  virtual void Reset() {};

 protected:
  int id_;

  int num_variables_;
  int num_reactors_;

  int num_species_;
  int num_steps_;
};

int ReactorGetTimeDerivative(const double reactor_time,
                      N_Vector state,
                      N_Vector derivative,
                      void* user_data);

#ifdef SUNDIALS2
int ReactorGetJacobianDense(long int N, double t, N_Vector y, N_Vector fy,
                     DlsMat Jac, void* user_data,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorGetJacobianDense(double t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                            void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#else
#error "Unsupported SUNDIALS version"
#endif

int ReactorGetJacobianDenseRaw(long int N, double t, N_Vector y, N_Vector fy,
                     double* Jac, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int ReactorJacobianSetup(double t, N_Vector y, N_Vector fy, void* user_data);

int ReactorJacobianFactor(realtype gamma, void* user_data);

#ifdef SUNDIALS2
int ReactorJacobianSetupAndFactor(double t, N_Vector y, N_Vector fy,
                                  booleantype jok, booleantype *jcurPtr,
                                  double gamma, void* user_data,
                                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorJacobianSetupAndFactor(double t, N_Vector y, N_Vector fy,
                                  booleantype jok, booleantype *jcurPtr,
                                  double gamma, void* user_data);
#else
#error "Unsupported SUNDIALS version"
#endif

#ifdef SUNDIALS2
int ReactorJacobianSolveCVODE(double t, N_Vector y, N_Vector fy,
                         N_Vector r, N_Vector z, realtype gamma,
                         realtype delta, int lr, void *user_data,
                         N_Vector tmp);
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorJacobianSolveCVODE(double t, N_Vector y, N_Vector fy,
                         N_Vector r, N_Vector z, realtype gamma,
                         realtype delta, int lr, void *user_data);
#else
#error "Unsupported SUNDIALS version"
#endif

int ReactorJacobianSolve(double t, N_Vector y, N_Vector fy,
                         N_Vector r, N_Vector z, void *user_data);

int ReactorComplexJacobianFactor(int k, double alpha, double beta, void* user_data);

int ReactorComplexJacobianSolve(int k, N_Vector ax, N_Vector bx, void* user_data);

int ReactorGetNumRootFunctions(void *user_data);

int ReactorRootFunction(double t, N_Vector y, double *root_function, void *user_data);
#endif
