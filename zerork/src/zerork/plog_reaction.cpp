#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "plog_reaction.h"

namespace zerork {

PLogReaction::PLogReaction(const int step_id,
                           const int num_pts,
                           const double pres[],
                           const double log_e_A[],
                           const double Tpow[],
                           const double Tact[],
                           const bool use_extrapolation)
{
  step_index_=step_id;
  num_pressure_points_=num_pts;
  use_extrapolation_=use_extrapolation;

  // initialize arrhenius parameter vectors to zero
  pressure_points_.assign(num_pts,0.0);
  log_e_pressure_points_.assign(num_pts,0.0);
  log_e_Afactor_.assign(num_pts,0.0);
  temperature_power_.assign(num_pts,0.0);
  activation_temperature_.assign(num_pts,0.0);

  // assign arrhenius parameter vectors
  for(int j=0; j<num_pts; ++j) {
    pressure_points_[j]=pres[j];
    log_e_Afactor_[j]     = log_e_A[j];
    temperature_power_[j] = Tpow[j];
    activation_temperature_[j] = Tact[j];
  }

  SortByPressure();
  for(int j=0; j<num_pts; ++j) {
    log_e_pressure_points_[j]=log(pressure_points_[j]);
  }

  max_pressure_ = pressure_points_[num_pts-1];
  min_pressure_ = pressure_points_[0];

  last_bounded_range_id_=0;
}
 
void PLogReaction::SortByPressure()
{
  int min_pressure_id;
  double min_pressure;
  double swap_value;

  for(int j=0; j<num_pressure_points_-1; ++j) {
    // find the minimum between j and n-2
    min_pressure_id = j;
    min_pressure    = pressure_points_[j];
    for(int k=j+1; k<num_pressure_points_; ++k) {
      if(pressure_points_[k] < min_pressure) {
        min_pressure_id = k;
        min_pressure    = pressure_points_[k];
      }
    }
    if(min_pressure_id != j) {
      // swap all parameters between min_pressure_id and j
      swap_value = pressure_points_[j];
      pressure_points_[j] = pressure_points_[min_pressure_id];
      pressure_points_[min_pressure_id] = swap_value;

      swap_value = log_e_Afactor_[j];
      log_e_Afactor_[j] = log_e_Afactor_[min_pressure_id];
      log_e_Afactor_[min_pressure_id] = swap_value;

      swap_value = temperature_power_[j];
      temperature_power_[j] = temperature_power_[min_pressure_id];
      temperature_power_[min_pressure_id] = swap_value;

      swap_value = activation_temperature_[j];
      activation_temperature_[j] = activation_temperature_[min_pressure_id];
      activation_temperature_[min_pressure_id] = swap_value;      
    }
  }
}
// GetPressureRangeIndex returns the vector index 'id' such that the
// PLOG data from pressure_points_[id] and pressure_points_[id+1] are to
// be used in GetRateConstantLogE(...).  If use_extrapolation_ is false and
// pressure is out of the range, then the 'id' returned is out of range to
// indicate special processing.  If use_extrapolation_ is true, 'id' is set
// to the first or last pressure range index. 
int PLogReaction::GetPressureRangeIndex(const double pressure)
{
  const int num_ranges=num_pressure_points_-1;

  if(min_pressure_ <= pressure && pressure < max_pressure_) {
    // check if the pressure is in the same range
    if(pressure_points_[last_bounded_range_id_] <= pressure &&
            pressure < pressure_points_[last_bounded_range_id_+1]) {
      return last_bounded_range_id_;
    }
    // otherwise scan all the ranges
    for(int j=0; j<num_ranges; ++j) {
      if(pressure_points_[j] <= pressure &&
         pressure < pressure_points_[j+1]) {
        return j;
      }
    }
    // this point is only reached in error
    printf("ERROR: In PLogReaction::GetPressureRangeIndex(...),\n");
    printf("       could not find pressure range index for p = %.18g [Pa]\n",
           pressure);
    printf("       in [%.18g, %.18g)\n",
           min_pressure_,
	   max_pressure_);
    exit(-1);

  } // end if(min <= pressure < max) 
  else if (pressure < min_pressure_) {
    return ((use_extrapolation_) ? 0 : -1);
  } else {
    return ((use_extrapolation_) ? num_ranges-1 : num_ranges);
  }
}


double PLogReaction::GetRateConstantLogE(const double inv_temperature,
                                         const double log_e_temperature,
                                         const double pressure,
				         const double log_e_pressure)
{
  double log_e_K,log_e_K_pt1, log_e_K_pt2 ;
  const int num_ranges = num_pressure_points_-1;
  int range_id = GetPressureRangeIndex(pressure);

  if(0 <= range_id && range_id < num_ranges) {

    // logarithmic interpolation (or extrapolation)
    log_e_K_pt1 = log_e_Afactor_[range_id] + 
      temperature_power_[range_id]*log_e_temperature -
      activation_temperature_[range_id]*inv_temperature;
    log_e_K_pt2 = log_e_Afactor_[range_id+1] + 
      temperature_power_[range_id+1]*log_e_temperature -
      activation_temperature_[range_id+1]*inv_temperature;

    log_e_K = log_e_K_pt1 + (log_e_K_pt2 - log_e_K_pt1)*
      (log_e_pressure - log_e_pressure_points_[range_id])/
      (log_e_pressure_points_[range_id+1]-log_e_pressure_points_[range_id]);

    last_bounded_range_id_ = range_id;
  } else if(range_id < 0) {
    // (pressure < pressure_min_ && use_extrapolation == false)
    log_e_K = log_e_Afactor_[0] + 
      temperature_power_[0]*log_e_temperature -
      activation_temperature_[0]*inv_temperature;
  } else {
    // (pressure >= pressure_max_ && use_extrapolation == false)
    log_e_K = log_e_Afactor_[num_ranges] + 
      temperature_power_[num_ranges]*log_e_temperature -
      activation_temperature_[num_ranges]*inv_temperature;
  } 
  return log_e_K;
}
double PLogReaction::GetRateConstantFromTP(const double temperature,
                                           const double pressure)
{
  return exp(GetRateConstantLogE(1.0/temperature,
                                 log(temperature),
                                 pressure,
                                 log(pressure)));
}


} // end of namespace zerork
