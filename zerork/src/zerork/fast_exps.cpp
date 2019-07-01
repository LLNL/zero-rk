#include <fast_exps.h>
#include <math.h> //for libc exp

#ifdef ZERORK_USE_AMD_LIBM
#include "amdlibm.h"
#endif

#ifdef ZERORK_USE_MKL
#include <mkl.h>
#endif

#if defined(ZERORK_USE_FMATH) || defined(ZERORK_USE_FMATH_NOCHECK)
#include <fmath.hpp>
#endif

namespace zerork {


void vecexp_libc(double *values, int n)
{
    int i;
    for (i = 0;i < n;++i) {
        values[i] = exp(values[i]);
    }
}

#ifdef ZERORK_USE_MKL
void vecexp_MKL_LA(double *values, int num)
{
  ::vmdExp(num,values,values,VML_LA);
}

void vecexp_MKL_HA(double *values, int num)
{
  ::vmdExp(num,values,values,VML_HA);
}

double exp_MKL_LA(double val)
{
  double out;
  ::vmdExp(1,&val,&out,VML_LA);
  return out;
}

double exp_MKL_HA(double val)
{
  double out;
  ::vmdExp(1,&val,&out,VML_HA);
  return out;
}
#endif


#ifdef ZERORK_USE_AMD_LIBM
void vecexp_AMD_LibM(double *values, int num)
{
  ::amd_vrda_exp(num,values,values);
}
#endif

#ifdef ZERORK_USE_FMATH
const double MIN_EXP_ARG = log(1.0e-300);
const double MAX_EXP_ARG = log(1.0e+300);
void vecexp_FMATH(double * x, int n) {
    for(int i = 0; i < n; i++) {
      x[i] = x[i] > MAX_EXP_ARG ? MAX_EXP_ARG : x[i];
      x[i] = x[i] < MIN_EXP_ARG ? MIN_EXP_ARG : x[i];
    }
    fmath::expd_v(x,n);
}
double exp_FMATH(double x) {
    x = x > MAX_EXP_ARG ? MAX_EXP_ARG : x;
    x = x > MIN_EXP_ARG ? MIN_EXP_ARG : x;
    return fmath::expd(x);
}
#endif

#ifdef ZERORK_USE_MKL
  #ifdef ZERORK_USE_MKL_HA
    void (*fast_vec_exp)(double *, int) = &(vecexp_MKL_HA);
    double (*fast_exp)(double) = &exp_MKL_HA;
  #else
    void (*fast_vec_exp)(double *, int) = &(vecexp_MKL_LA);
    double (*fast_exp)(double) = &exp_MKL_LA;
  #endif
#else
  #ifdef ZERORK_USE_AMD_LIBM
    void (*fast_vec_exp)(double *, int) = &(vecexp_AMD_LibM);
    double (*fast_exp)(double) = &(exp); //AMD_LibM overrides ::exp
  #else
    #ifdef ZERORK_USE_FMATH
      void (*fast_vec_exp)(double *, int) = &(vecexp_FMATH);
//      double (*fast_exp)(double) = &(exp_FMATH);
      double (*fast_exp)(double) = &(::exp); //TODO:fmath::expd gives bad values
    #else
      #ifdef ZERORK_USE_FMATH_NOCHECK
        // No test to ensure in range
        void (*fast_vec_exp)(double *, int) = &(fmath::expd_v);
//        double (*fast_exp)(double) = &(fmath::expd);
        double (*fast_exp)(double) = &(::exp); //TODO:fmath::expd gives bad values
      #else
   // USE LIBC
        void (*fast_vec_exp)(double *, int) = &(vecexp_libc);
        double (*fast_exp)(double) = &(::exp);
      #endif
    #endif
  #endif
#endif


} // namespace zerork
