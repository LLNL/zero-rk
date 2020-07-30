#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#include "gsa_stats.h"

GSAStats::GSAStats(const int dim)
{
  num_dimensions_ = dim;
  num_samples_ = 0;

  min_value_.assign(dim, 1.0e300);
  max_value_.assign(dim,-1.0e300);
  log_mean_sum_.assign(dim, 0.0);
  log_variance_sum_.assign(dim, 0.0);
}

void GSAStats::AddNextSample(const double sample_data[])
{
  const int const_num_dimensions = num_dimensions_;

  for(int j=0; j<const_num_dimensions; ++j) {
    if(sample_data[j] < min_value_[j]) {
      min_value_[j] = sample_data[j];
    }
    if(sample_data[j] > max_value_[j]) {
      max_value_[j] = sample_data[j];
    }
    double log_sample = log(sample_data[j]);
    log_mean_sum_[j] += log_sample;
    log_variance_sum_[j] += log_sample*log_sample;
  }
  ++num_samples_;
}


bool GSAStats::ValidIndex(const int dim) const 
{
  if(0 <= dim && dim <num_dimensions_) {
    return true;
  }
  printf("WARNING: Invalid index %d for array of length %d\n",
         dim,num_dimensions_);
  fflush(stdout);
  return false;
}

 
double GSAStats::MinValue(const int dim) const
{
  if(ValidIndex(dim)) {
    return min_value_[dim];
  }
  return 0.0;
}
double GSAStats::MaxValue(const int dim) const
{
  if(ValidIndex(dim)) {
    return max_value_[dim];
  }
  return 0.0;
}
double GSAStats::LogMean(const int dim) const
{
  if(ValidIndex(dim)) {
    if(num_samples_ <= 0) {
      printf("WARNING: In GSAStats::Mean(...)\n");
      printf("         num_samples = %d\n",num_samples_);
      return 0.0; 
    }
    return log_mean_sum_[dim]/(double)num_samples_;
  }
  return 0.0;
}
double GSAStats::LogVariance(const int dim) const
{
  double variance_dim, mean_dim;
  if(ValidIndex(dim)) {
    if(num_samples_ <= 1) {
      printf("WARNING: In GSAStats::Variance(...)\n");
      printf("         num_samples = %d\n",num_samples_); 
      return 0.0;
    }
    mean_dim = LogMean(dim);
    variance_dim = log_variance_sum_[dim]/(double)num_samples_ - 
      mean_dim*mean_dim;

    return variance_dim;
  }
  return 0.0;
}



