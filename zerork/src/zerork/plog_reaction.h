#ifndef PLOG_REACTION_H_
#define PLOG_REACTION_H_

#include <vector>

namespace zerork {


class PLogReaction {

 public:
  PLogReaction(const int step_id,
               const int num_pts,
               const double pres[],
               const double log_e_A[],
               const double Tpow[],
               const double Tact[],
               const bool use_extrapolation);
  //~PLogReaction();
  int num_pressure_points() const {return num_pressure_points_;}
  bool use_extrapolation() const {return use_extrapolation_;}
  double max_pressure() const {return max_pressure_;}
  double min_pressure() const {return min_pressure_;}
  int step_index() const {return step_index_;}
  int GetPressureRangeIndex(const double pressure);
  double GetRateConstantLogE(const double inv_temperature,
                             const double log_e_temperature,
                             const double pressure,
                             const double log_e_pressure);
  double GetRateConstantFromTP(const double temperature,
                               const double pressure);
  
 private:
  void SortByPressure();

  double max_pressure_;
  double min_pressure_;
  int num_pressure_points_;
  int last_bounded_range_id_;
  int step_index_;
  bool use_extrapolation_;

  std::vector<double> pressure_points_;
  std::vector<double> log_e_pressure_points_;
  std::vector<double> log_e_Afactor_;
  std::vector<double> temperature_power_;
  std::vector<double> activation_temperature_;
};

} // end of namespace zerork
#endif
