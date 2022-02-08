
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
                     DlsMat Jac, void* user_data,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->GetJacobianDense(N, t, y, fy, Jac, tmp1, tmp2, tmp3);
}
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorGetJacobianDense(double t, N_Vector y, N_Vector fy,
                            SUNMatrix Jac, void* user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->GetJacobianDense(t, y, fy, Jac, tmp1, tmp2, tmp3);
}
#else
#error "Unsupported SUNDIALS version"
#endif

int ReactorJacobianSetup(double t, N_Vector y, N_Vector fy, void* user_data,
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->JacobianSetup(t, y, fy, tmp1, tmp2, tmp3);
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
  return reactor->JacobianSolve(t, y, fy, r, z, tmp);
}
#elif defined SUNDIALS3 || defined SUNDIALS4
int ReactorJacobianSolveCVODE(double t, N_Vector y, N_Vector fy,
                         N_Vector r, N_Vector z, realtype gamma,
                         realtype delta, int lr, void *user_data)
{
  N_Vector tmp;
  tmp = N_VClone(y);
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  int flag = reactor->JacobianSolve(t, y, fy, r, z, tmp);
  N_VDestroy(tmp);
  return flag;
}
#else
#error "Unsupported SUNDIALS version"
#endif

int ReactorJacobianSolve(double t, N_Vector y, N_Vector fy,
                         N_Vector r, N_Vector z, void *user_data,
                         N_Vector tmp)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  return reactor->JacobianSolve(t, y, fy, r, z, tmp);
}

#ifdef SUNDIALS2
int ReactorJacobianSetupAndFactor(double t, N_Vector y, N_Vector fy,
                                  booleantype jok, booleantype *jcurPtr,
                                  double gamma, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  if(!jok) {
    (*jcurPtr)=TRUE;
    reactor->JacobianSetup(t,y,fy,tmp1,tmp2,tmp3);
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
  N_Vector tmp1, tmp2, tmp3;
  tmp1 = N_VClone(y);
  tmp2 = N_VClone(y);
  tmp3 = N_VClone(y);
  ReactorBase* reactor = static_cast<ReactorBase*>(user_data);
  if(!jok) {
    (*jcurPtr)=SUNTRUE;
    reactor->JacobianSetup(t,y,fy,tmp1,tmp2,tmp3);
  } else {
    (*jcurPtr)=SUNFALSE;
  }
  int flag = reactor->JacobianFactor(gamma);
  N_VDestroy(tmp1);
  N_VDestroy(tmp2);
  N_VDestroy(tmp3);
  return flag;
}
#else
#error "Unsupported SUNDIALS version"
#endif

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

