#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

#include <vector>

namespace zerork {
namespace utilities {
class Distribution
{
 public:
  Distribution(const int n,
               const double x_min, 
               const double x_max,
               const bool use_log);
  Distribution(); // default distribution is linear, [0,1) with one bin
  // copy constructor
  Distribution(const Distribution &dist);
  // basic operations
  Distribution& operator=(const Distribution &rdist);
  Distribution& operator+=(const Distribution &rdist);
  Distribution  operator+(const Distribution &rdist);

  // basic accessors
  bool use_log_spacing() const {return use_log_spacing_;}
  int num_bins() const {return num_bins_;}
  double total_weight() const {return total_weight_;}
  double under_range_weight() const {return under_range_weight_;}
  double over_range_weight() const {return over_range_weight_;}
  double min_range() const {return min_range_;}
  double max_range() const {return max_range_;}
  // special accessors
  int getBinIndex(const double val) const;
  double getBinWeight(const int id) const;
  const double *getBinWeightData() const {return bin_weight_.data();}
  double getBinFraction(const int id) const
  {return getBinWeight(id)/total_weight_;}
  double getBinMin(const int id) const;
  double getBinMax(const int id) const;
  const double *getBinBoundaryData() const {return bin_boundary_.data();}

  int addValue(const double val, const double wt);
  bool isSameBinDistributionAs(const Distribution &rdist);
  void zeroBinWeights();

 private:
  bool isBinIndexInRange(const int id) const;
  bool isValueInRange(const double val) const;
  int LinearBinMap(const double val) const;
  int LogBinMap(const double val) const;

  bool use_log_spacing_;
  int num_bins_;
  double min_range_;
  double max_range_;
  double total_weight_;
  double under_range_weight_;
  double over_range_weight_;
  double log_min_range_;
  double log_delta_range_;
  std::vector<double> bin_weight_;
  std::vector<double> bin_boundary_;

};
} // end of namespace utilities
} // end of namespace zerork
#endif
