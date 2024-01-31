#include "fast_exps.h"

#include "zerork_conf.h"

#ifndef ZERORK_EXP_LIBC
  #ifdef ZERORK_HAVE_MKL
    #define ZERORK_EXP_MKL
  #else
    #ifdef ZERORK_HAVE_IBM_MASS
      #define ZERORK_EXP_MASS
    #else
      #define ZERORK_EXP_FMATH
    #endif
  #endif
#endif

#include <math.h> //for libc exp

#ifdef ZERORK_EXP_AMD_LIBM
#include "amdlibm.h"
#endif

#ifdef ZERORK_EXP_MKL
#include "mkl.h"
#endif

#if defined(ZERORK_EXP_FMATH) || defined(ZERORK_EXP_FMATH_NOCHECK)
#include <fmath.hpp>
#include "zerork/constants.h"
#endif

#ifdef ZERORK_EXP_MASS
#include <massv.h>
#endif


namespace zerork {


void vecexp_libc(double *values, int n)
{
    int i;
    for (i = 0;i < n;++i) {
        values[i] = exp(values[i]);
    }
}

#ifdef ZERORK_EXP_MKL
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


#ifdef ZERORK_EXP_AMD_LIBM
void vecexp_AMD_LibM(double *values, int num)
{
  ::amd_vrda_exp(num,values,values);
}
#endif

#ifdef ZERORK_EXP_FMATH
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

#ifdef ZERORK_EXP_MASS
void vecexp_MASS(double *values, int num)
{
  ::vexp(values, values, &num);
}

double exp_MASS(double val)
{
  double out;
  int num = 1;
  ::vexp(&out,&val,&num);
  return out;
}
#endif

#ifdef ZERORK_EXP_MKL
  #ifdef ZERORK_EXP_MKL_HA
    const char* expType="intel_mkl_ha";
    void (*fast_vec_exp)(double *, int) = &(vecexp_MKL_HA);
    double (*fast_exp)(double) = &exp_MKL_HA;
  #else
    const char* expType="intel_mkl_la";
    void (*fast_vec_exp)(double *, int) = &(vecexp_MKL_LA);
    double (*fast_exp)(double) = &exp_MKL_LA;
  #endif
#else
  #ifdef ZERORK_EXP_AMD_LIBM
    const char* expType="amd_libm";
    void (*fast_vec_exp)(double *, int) = &(vecexp_AMD_LibM);
    double (*fast_exp)(double) = &(exp); //AMD_LibM overrides ::exp
  #else
    #ifdef ZERORK_EXP_MASS
      const char* expType="ibm_mass";
      void (*fast_vec_exp)(double *, int) = &(vecexp_MASS);
      double (*fast_exp)(double) = &exp_MASS;
    #else
      #ifdef ZERORK_EXP_FMATH
        const char* expType="fmath";
        void (*fast_vec_exp)(double *, int) = &(vecexp_FMATH);
//      double (*fast_exp)(double) = &(exp_FMATH);
        double (*fast_exp)(double) = &(::exp); //TODO:fmath::expd gives bad values
      #else
        #ifdef ZERORK_EXP_FMATH_NOCHECK
          // No test to ensure in range
          const char* expType="fmath_nocheck";
          void (*fast_vec_exp)(double *, int) = &(fmath::expd_v);
//        double (*fast_exp)(double) = &(fmath::expd);
          double (*fast_exp)(double) = &(::exp); //TODO:fmath::expd gives bad values
        #else
     // ZERORK_EXP LIBC
          const char* expType="libc";
          void (*fast_vec_exp)(double *, int) = &(vecexp_libc);
          double (*fast_exp)(double) = &(::exp);
        #endif
      #endif
    #endif
  #endif
#endif


} // namespace zerork
