
#include "reactor_base.h"

int ReactorGetTimeDerivative(const double reactor_time,
                      N_Vector state,
                      N_Vector derivative,
                      void* user_data)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->GetTimeDerivative(reactor_time, state, derivative);
}

#ifdef SUNDIALS2
int ReactorGetJacobianDense(long int N, double t, N_Vector y, N_Vector fy,
                     DlsMat Jac, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->GetJacobianDense(N, t, y, fy, Jac);
}
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorGetJacobianDense(double t, N_Vector y, N_Vector fy,
                            SUNMatrix Jac, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->GetJacobianDense(t, y, fy, Jac);
}
#else
#error "Unsupported SUNDIALS version"
#endif

int ReactorGetJacobianDenseRaw(long int N, double t, N_Vector y, N_Vector fy,
                               double* Jac, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->GetJacobianDenseRaw(N, t, y, fy, Jac);
}

int ReactorJacobianSetup(double t, N_Vector y, N_Vector fy, void* user_data)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->JacobianSetup(t, y, fy);
}

int ReactorJacobianFactor(double gamma, void* user_data)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->JacobianFactor(gamma);
}

#ifdef SUNDIALS2
int ReactorJacobianSolveCVODE(double t, N_Vector y, N_Vector fy,
                         N_Vector r, N_Vector z, realtype gamma,
                         realtype delta, int lr, void *user_data,
                         N_Vector tmp)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->JacobianSolve(t, y, fy, r, z);
}
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorJacobianSolveCVODE(double t, N_Vector y, N_Vector fy,
                         N_Vector r, N_Vector z, realtype gamma,
                         realtype delta, int lr, void *user_data)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  int flag = reactor->JacobianSolve(t, y, fy, r, z);
  return flag;
}
#else
#error "Unsupported SUNDIALS version"
#endif

int ReactorJacobianSolve(double t, N_Vector y, N_Vector fy,
                         N_Vector r, N_Vector z, void *user_data)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->JacobianSolve(t, y, fy, r, z);
}

#ifdef SUNDIALS2
int ReactorJacobianSetupAndFactor(double t, N_Vector y, N_Vector fy,
                                  booleantype jok, booleantype *jcurPtr,
                                  double gamma, void* user_data,
                                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  if(!jok) {
    (*jcurPtr)=TRUE;
    reactor->JacobianSetup(t,y,fy);
  } else {
    (*jcurPtr)=FALSE;
  }
  int flag = reactor->JacobianFactor(gamma);
  return flag;
}
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorJacobianSetupAndFactor(double t, N_Vector y, N_Vector fy,
                                  booleantype jok, booleantype *jcurPtr,
                                  double gamma, void* user_data)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  if(!jok) {
    (*jcurPtr)=SUNTRUE;
    reactor->JacobianSetup(t,y,fy);
  } else {
    (*jcurPtr)=SUNFALSE;
  }
  int flag = reactor->JacobianFactor(gamma);
  return flag;
}
#else
#error "Unsupported SUNDIALS version"
#endif

int ReactorComplexJacobianFactor(int k, double alpha, double beta, void* user_data)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->ComplexJacobianFactor(k, alpha, beta);
}

int ReactorComplexJacobianSolve(int k, N_Vector ax, N_Vector bx, void *user_data)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->ComplexJacobianSolve(k, ax, bx);
}

int ReactorGetNumRootFunctions(void *user_data)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->GetNumRootFunctions();
}

int ReactorRootFunction(double t, N_Vector y, double *root_function,
                        void *user_data)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->RootFunction(t, y, root_function);
}

int ReactorBase::GetID() {
 return id_;
}

void ReactorBase::SetID(int id) {
  id_ = id;
}

