#ifndef PLOG_REACTION_H_
#define PLOG_REACTION_H_

#include <vector>

namespace zerork {

class rate_const_cuda; //forward declaration for friendship

class PLogReaction {
  friend class zerork::rate_const_cuda;//needs direct access to tables

 public:
  PLogReaction(const int rxn_id,
               const int step_id,
               const int num_pts,
               const double pres[],
               const double Afactor[],
               const double Tpow[],
               const double Tact[],
               const bool use_extrapolation);

  //~PLogReaction();
  int total_arrhenius_lines() const {return total_arrhenius_lines_;}
  int num_pressure_points() const {return num_pressure_points_;}
  int num_arrhenius_lines(const int pressure_id) const {return num_lines_at_pressure_[pressure_id];}
  bool use_extrapolation() const {return use_extrapolation_;}
  double max_pressure() const {return max_pressure_;}
  double min_pressure() const {return min_pressure_;}
  int reaction_index() const {return reaction_index_;}
  int step_index() const {return step_index_;}
  int GetPressureRangeIndex(const double pressure);
  double GetRateCoefficientFromTP(const double tempreature,
                               const double inv_temperature,
                               const double log_e_temperature,
                               const double pressure,
                               const double log_e_pressure);
  double GetRateCoefficientFromTP(const double temperature,
                               const double pressure);
  
 private:
  void SortByPressure();
  void CountDuplicatePressures();
  double GetRateCoefficientAtPressure(const int pressure_id,
                                      const double temperature,
                                      const double inv_temperature,
                                      const double log_e_temperature); 


  double max_pressure_;
  double min_pressure_;
  int num_pressure_points_;
  int total_arrhenius_lines_;
   
  int last_bounded_range_id_;
  int reaction_index_;
  int step_index_;
  bool use_extrapolation_;

  std::vector<int> num_lines_at_pressure_;
  std::vector<int> start_line_at_pressure_;
  std::vector<double> pressure_points_;
  std::vector<double> log_e_pressure_points_;
  std::vector<double> log_e_afactor_;
  std::vector<double> sign_afactor_;
  std::vector<double> temperature_power_;
  std::vector<double> activation_temperature_;
};

} // end of namespace zerork
#endif
