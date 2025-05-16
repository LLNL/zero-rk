#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "plog_reaction.h"

namespace zerork {

static const double PRESSURE_RTOL = 1.0e-8;

PLogReaction::PLogReaction(const int rxn_id,
                           const int step_id,
                           const int num_pts,
                           const double pres[],
                           const double Afactor[],
                           const double Tpow[],
                           const double Tact[],
                           const bool use_extrapolation)
{
  reaction_index_=rxn_id;
  step_index_=step_id;
  total_arrhenius_lines_=num_pts;
  use_extrapolation_=use_extrapolation;

  // initialize arrhenius parameter vectors to zero
  pressure_points_.assign(num_pts,0.0);
  log_e_afactor_.assign(num_pts,0.0);
  sign_afactor_.assign(num_pts, 0.0);
  temperature_power_.assign(num_pts,0.0);
  activation_temperature_.assign(num_pts,0.0);

  // assign arrhenius parameter vectors
  for(int j=0; j<num_pts; ++j) {
    pressure_points_[j]        = pres[j];
    log_e_afactor_[j]          = log(fabs(Afactor[j]));
    sign_afactor_[j]           = ((Afactor[j] < 0.0) ? -1.0 : 1.0);
    temperature_power_[j]      = Tpow[j];
    activation_temperature_[j] = Tact[j];
  }

  SortByPressure();
  CountDuplicatePressures();

  log_e_pressure_points_.clear();
  for(int j=0; j<num_pressure_points_; ++j) {
    log_e_pressure_points_.push_back(log(pressure_points_[j]));
  }

  max_pressure_ = pressure_points_[num_pressure_points_-1];
  min_pressure_ = pressure_points_[0];

  last_bounded_range_id_=0;
}
 
void PLogReaction::SortByPressure()
{
  int min_pressure_id;
  double min_pressure;
  double swap_value;

  for(int j=0; j<total_arrhenius_lines_-1; ++j) {
    // find the minimum between j and n-2
    min_pressure_id = j;
    min_pressure    = pressure_points_[j];
    for(int k=j+1; k<total_arrhenius_lines_; ++k) {
      if(pressure_points_[k] < min_pressure) {
        min_pressure_id = k;
        min_pressure    = pressure_points_[k];
      }
    }
    if(min_pressure_id != j) {
      // swap all parameters between min_pressure_id and j
      swap_value                        = pressure_points_[j];
      pressure_points_[j]               = pressure_points_[min_pressure_id];
      pressure_points_[min_pressure_id] = swap_value;

      swap_value                      = log_e_afactor_[j];
      log_e_afactor_[j]               = log_e_afactor_[min_pressure_id];
      log_e_afactor_[min_pressure_id] = swap_value;

      swap_value                     = sign_afactor_[j];
      sign_afactor_[j]               = sign_afactor_[min_pressure_id];
      sign_afactor_[min_pressure_id] = swap_value;

      swap_value = temperature_power_[j];
      temperature_power_[j] = temperature_power_[min_pressure_id];
      temperature_power_[min_pressure_id] = swap_value;

      swap_value = activation_temperature_[j];
      activation_temperature_[j] = activation_temperature_[min_pressure_id];
      activation_temperature_[min_pressure_id] = swap_value;      
    }
  }
}
// requires the pressures to be sorted, will reshrink the pressure_points_
// vector to the number of distinct pressure points
void PLogReaction::CountDuplicatePressures()
{
  double last_pressure;
  std::vector<double> distinct_pressures;
  int last_id;

  distinct_pressures.clear();
  num_lines_at_pressure_.clear();
  start_line_at_pressure_.clear();

  if(pressure_points_.size() > 0) {

    distinct_pressures.push_back(pressure_points_[0]);
    last_pressure = pressure_points_[0];

    num_lines_at_pressure_.push_back(1);
    last_id = 0;
    start_line_at_pressure_.push_back(0);

    for(size_t j=1; j<pressure_points_.size(); ++j) {

      if(fabs(pressure_points_[j]-last_pressure) > 
         fabs(last_pressure*PRESSURE_RTOL)) {
        // new pressure
        distinct_pressures.push_back(pressure_points_[j]);
        last_pressure = pressure_points_[j];
        num_lines_at_pressure_.push_back(1);
        ++last_id;
        start_line_at_pressure_.push_back(j); 

      } else {
        // duplicate pressure
        ++num_lines_at_pressure_[last_id];
      }

    }
    pressure_points_ = distinct_pressures;
    num_pressure_points_ = (int)pressure_points_.size();

    // TODO: add sanity check on vector sizes and ordering

  }
}

// GetPressureRangeIndex returns the vector index 'id' such that the
// PLOG data from pressure_points_[id] and pressure_points_[id+1] are to
// be used in GetRateCoefficientFromTP(...).  If use_extrapolation_ is false 
// and pressure is out of the range, then the 'id' returned is out of range to
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
    assert(("Couldn't find pressure range index", false));
    return -2; // avoid compiler warning

  } // end if(min <= pressure < max) 
  else if (pressure < min_pressure_) {
    return ((use_extrapolation_) ? 0 : -1);
  } else {
    return ((use_extrapolation_) ? num_ranges-1 : num_ranges);
  }
}

double PLogReaction::GetRateCoefficientAtPressure(const int pressure_id,
                                               const double temperature,
                                               const double inv_temperature,
                                               const double log_e_temperature)
{
  double rate_coefficient = 0.0;
 
  if(pressure_id < 0 || pressure_id >= (int)pressure_points_.size()) {
    return rate_coefficient;
  }

  const int num_lines = num_lines_at_pressure_[pressure_id];
  int line_id = start_line_at_pressure_[pressure_id];

  for(int j=0; j<num_lines; ++j) {
    rate_coefficient += sign_afactor_[line_id+j]*
      exp(log_e_afactor_[line_id+j] +
          log_e_temperature*temperature_power_[line_id+j] -
          inv_temperature*activation_temperature_[line_id+j]);
  }
  return rate_coefficient;
} 

double PLogReaction::GetRateCoefficientFromTP(const double temperature,
                                           const double inv_temperature,
                                           const double log_e_temperature,
                                           const double pressure,
				           const double log_e_pressure)
{
  double rate_coefficient_p1, rate_coefficient_p2;
  const int num_ranges = num_pressure_points_-1;
  double interpolation_exponent = 0.0;
  double rate_coefficient = 0.0;
  static bool neg_pressure_warned = false;
  static bool neg_coeff_warned = false;

  if(temperature <= 0.0 || pressure <= 0.0) {
    if(!neg_pressure_warned) {
      neg_pressure_warned=true;
      fprintf(stderr,"# WARNING: In PLogReaction::GetRateCoefficientFromTP(...),\n");
      fprintf(stderr,"#          can not use log interpolation on a negative pressure\n");
      fprintf(stderr,"#          or temperature:\n");
      fprintf(stderr,"#              p = %.18g [Pa]\n", pressure);
      fprintf(stderr,"#              T = %.18g [K]\n", temperature);
      fprintf(stderr,"#          Reaction %d\n", reaction_index_);
      fprintf(stderr,"#          Returning a rate coefficient of zero.\n");
    }
    return 0.0;
  }

  if(num_ranges == 0) {
    // single pressure
    return GetRateCoefficientAtPressure(0,  // pressure point index
                                        temperature,
                                        inv_temperature,
                                        log_e_temperature);
  }

  if(pressure <= min_pressure_ && use_extrapolation_ == false) {

    return GetRateCoefficientAtPressure(0,  // pressure point index
                                        temperature,
                                        inv_temperature,
                                        log_e_temperature);
    
  } else if(pressure >= max_pressure_ && use_extrapolation_ == false) {

    return GetRateCoefficientAtPressure(num_pressure_points_-1, 
                                        temperature,
                                        inv_temperature,
                                        log_e_temperature);
    
  }
  
  // compute the rate coefficient using logarithmic interploation
  int range_id = GetPressureRangeIndex(pressure);

  if(0 <= range_id && range_id < num_ranges) {

    last_bounded_range_id_ = range_id;

    // logarithmic interpolation (or extrapolation)
    rate_coefficient_p1 = GetRateCoefficientAtPressure(range_id, 
                                                       temperature,
                                                       inv_temperature,
                                                       log_e_temperature);
    if(rate_coefficient_p1 <= 0.0) {
      if(!neg_coeff_warned) {
        neg_coeff_warned = true;
        fprintf(stderr,"# ERROR:   In PLogReaction::GetRateCoefficientFromTP(...),\n");
        fprintf(stderr,"#          can not use log interpolation on a negative rate coefficient:\n");
        fprintf(stderr,"#              K(p,T) = %.18g\n",
               rate_coefficient_p1);
        fprintf(stderr,"#                  p  = %.18g [Pa]\n",
               pressure_points_[range_id]);
        fprintf(stderr,"#                  T  = %.18g [K]\n",
               temperature);
        fprintf(stderr,"#          Reaction %d\n", reaction_index_);
        fprintf(stderr,"#          Returning a rate coefficient of zero.\n");
      }
      return 0.0;
    }
    rate_coefficient_p2 = GetRateCoefficientAtPressure(range_id+1, 
                                                       temperature,
                                                       inv_temperature,
                                                       log_e_temperature);
    if(rate_coefficient_p2 <= 0.0) {
      if(!neg_coeff_warned) {
        neg_coeff_warned = true;
        fprintf(stderr,"# ERROR:   In PLogReaction::GetRateCoefficientFromTP(...),\n");
        fprintf(stderr,"#          can not use log interpolation on a negative rate coefficient:\n");
        fprintf(stderr,"#              K(p,T) = %.18g\n",
               rate_coefficient_p2);
        fprintf(stderr,"#                  p  = %.18g [Pa]\n",
               pressure_points_[range_id+1]);
        fprintf(stderr,"#                  T  = %.18g [K]\n",
               temperature);
        fprintf(stderr,"#          Reaction %d\n", reaction_index_);
        fprintf(stderr,"#          Returning a rate coefficient of zero.\n");
      }
      return 0.0;
    }

    interpolation_exponent = 
      log(rate_coefficient_p2/rate_coefficient_p1)/
      (log_e_pressure_points_[range_id+1]-log_e_pressure_points_[range_id]);

    rate_coefficient = 
      rate_coefficient_p1*pow(pressure/pressure_points_[range_id],
                              interpolation_exponent);

  } else if(range_id < 0) {
    // (pressure < pressure_min_ && use_extrapolation == false)
    // TODO: the code should not be here, clean up
    rate_coefficient = GetRateCoefficientAtPressure(0, // pressure point index
                                                    temperature,
                                                    inv_temperature,
                                                    log_e_temperature);
  } else {
    // (pressure >= pressure_max_ && use_extrapolation == false)
    // TODO: the code should not be here, clean up
    rate_coefficient = GetRateCoefficientAtPressure(num_pressure_points_-1, 
                                                    temperature,
                                                    inv_temperature,
                                                    log_e_temperature);
  } 
  return rate_coefficient;
}
double PLogReaction::GetRateCoefficientFromTP(const double temperature,
                                              const double pressure)
{
  if(temperature <= 0.0 || pressure <= 0.0) {
    fprintf(stderr,"# WARNING: In PLogReaction::GetRateCoefficientFromTP(...),\n");
    fprintf(stderr,"#          can not use log interpolation on a negative pressure\n");
    fprintf(stderr,"#          or temperature:\n");
    fprintf(stderr,"#              p = %.18g [Pa]\n", pressure);
    fprintf(stderr,"#              T = %.18g [K]\n", temperature);
    fprintf(stderr,"#          Reaction %d\n", reaction_index_);
    fprintf(stderr,"#          Returning a rate coefficient of zero.\n");
    return 0.0;
  }

  return GetRateCoefficientFromTP(temperature,
                                  1.0/temperature,
                                  log(temperature),
                                  pressure,
                                  log(pressure));
}


} // end of namespace zerork
