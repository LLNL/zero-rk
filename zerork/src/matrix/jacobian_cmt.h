#ifndef JACOBIAN_CMT_H_
#define JACOBIAN_CMT_H_

#include "zerork/mechanism.h"
#include "exploded_rop_derivative.h"

namespace zerork {

int ConstVolumeSparseJacobian_CMT(mechanism *mechp,
                                  ExplodedROPDerivative *exploded_deriv,
                                  const double conc_orig[],
                                  const double mix_orig,
                                  const double temp_orig,
                                  const double min_conc,
                                  const double state_deriv[],
                                  const int num_non_zeros,
                                  const int column_sum[],
                                  const int row_id[], 
                                  double workspace[],
                                  double jacobian[]);



} // namespace zerork
#endif
