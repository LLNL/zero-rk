#ifndef MULTIZONE_H_
#define MULTIZONE_H_


#include "zerork/mechanism.h"

void multizone_make_zones(const zerork::mechanism &mech,
                          int nreactors, const double *T, const double *P,
                          const double *dpdt, const double *volume,
                          const double *massfracs, const double *cost,
                          const double *gpu,
                          int &nzones,
                          std::vector<double> &temp_zones,
                          std::vector<double> &press_zones,
                          std::vector<double> &dpdt_zones,
                          std::vector<double> &massfracs_zones,
                          std::vector<double> &cost_zones,
                          std::vector<double> &gpu_flag_zones);


void multizone_remap(const zerork::mechanism &mech,
                     int nzones,
                     const std::vector<double> &massfracs_zones0,
                     const std::vector<double> &massfracs_zones,
                     int nreactors,
                     double* massfracs,
                     const std::vector<double> &cost_zones,
                     double* cost,
                     const std::vector<double> &gpu_flag_zones,
                     double* gpu);
#endif
