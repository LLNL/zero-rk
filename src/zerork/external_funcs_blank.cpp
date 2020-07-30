#include "external_funcs.h"
#include <cstdio>
#include <cstdlib>

namespace zerork {

#ifdef EXIT_THROWS_EXCEPTION
  // create a local function to overide the system exit and throw an exception
  // with the status integer.
  static void exit(int status) {throw status;}
#endif // EXIT_THROWS_EXCEPTION

void external_func_check(int nsp, int nstep)
{
    std::printf("Attempt to use blank external func.  Quitting.\n");
#ifndef EXIT_THROWS_EXCEPTION
    std::exit(1);
#else
    exit(1);
#endif
}

void external_func_rates(const double C[], double phi[], double cdot[], double ddot[])
{}

void external_func_arrh(const double Tcurrent, double expWorkArray[], double Kwork[],
                            const int nDistinctArrhenius,
                            const double distinctArrheniusLogAfact[],
                            const double distinctArrheniusTpow[],
                            const double distinctArrheniusTact[])
{}

void external_func_keq(const int nFromKeqStep, const double Gibbs_RT[], double expWorkArray[], double Kwork[], const double log_e_PatmInvRuT)
{}

} // namespace zerork

