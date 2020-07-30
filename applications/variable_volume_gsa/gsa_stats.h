#ifndef GSA_STATS_H_
#define GSA_STATS_H_

#include <vector>

class GSAStats {
 public:
  GSAStats(const int dim);
  ~GSAStats() {0;}

  void AddNextSample(const double sample_data[]);

  double MinValue(const int dim) const;
  double MaxValue(const int dim) const;
  double LogMean(const int dim) const;
  double LogVariance(const int dim) const;

  // basic accessors
  int num_dimensions() const {return num_dimensions_;}
  int num_samples() const {return num_samples_;}

 private:
  bool ValidIndex(const int dim) const; 

  int num_dimensions_;
  int num_samples_;

  std::vector<double> min_value_;
  std::vector<double> max_value_;
  std::vector<double> log_mean_sum_;
  std::vector<double> log_variance_sum_;
  
};

#endif
