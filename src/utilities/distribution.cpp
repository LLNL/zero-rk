#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <limits>

#include "distribution.h"

namespace zerork {
namespace utilities {
Distribution::Distribution(const int n,
                           const double x_min, 
                           const double x_max,
                           const bool use_log)
{
  num_bins_        = n;
  min_range_       = x_min;
  max_range_       = x_max;
  use_log_spacing_ = use_log;
  //check to make sure meaning values are assigned in the constructor
  if(n < 1) {
    printf("ERROR: In Distribution::Distribution(...),\n");
    printf("       distribution must contain at least one bin\n");
    printf("       instead of n = %d\n",n);
    exit(-1);
  }
  if(min_range_ >= max_range_) {
    printf("ERROR: In Distribution::Distribution(...),\n");
    printf("       distribution must span a non-zero, increasing range.\n");
    printf("       instead of range = [%.18g,%.18g)\n",
           min_range_,max_range_);
    exit(-1);
  }

  log_delta_range_ = 0.0;
  log_min_range_   = 0.0;
  if(use_log_spacing_) {
    log_delta_range_  = log(max_range_/min_range_); // log(max) - log(min)
    log_min_range_    = log(min_range_);
  }  // should only be used if use_log_spacing_ = true

  // initialize weight bins and bin boundaries
  total_weight_       = 0.0;
  under_range_weight_ = 0.0;
  over_range_weight_  = 0.0;
  bin_weight_.assign(num_bins_,0.0);
  bin_boundary_.assign(num_bins_+1,0.0);
 
  // Assign the boundaries of the bin. Note that the boundaries of bin[i] are
  // stored in bin_boundary_[i] (min) and bin_boundary_[i+1] (max). 
  bin_boundary_[0] = min_range_;
  bin_boundary_[num_bins_] = max_range_;
  if(use_log_spacing_) {
    if(x_min <= 0.0) {
      printf("ERROR: In Distribution::Distribution(...),\n");
      printf("       logarithmic spacing requires a positive range min,\n");
      printf("       range min = %.18g\n",min_range_);
      exit(-1);
    }
    double delta_log;
    for(int j=1; j<num_bins_; ++j) {
      delta_log = static_cast<double>(j)/static_cast<double>(num_bins_)*
        log_delta_range_;
      bin_boundary_[j] = exp(log_min_range_ + delta_log);
    }
  } else { // use linear spacing
    double delta_linear;
    for(int j=1; j<num_bins_; ++j) {
      delta_linear = static_cast<double>(j)/static_cast<double>(num_bins_)*
        (max_range_-min_range_);
      bin_boundary_[j] = min_range_+delta_linear;
    }
  }
}

Distribution::Distribution()
{
  use_log_spacing_    = false;
  num_bins_           = 1;
  min_range_          = 0.0;
  max_range_          = 1.0;
  total_weight_       = 0.0;
  under_range_weight_ = 0.0;
  over_range_weight_  = 0.0;
  log_min_range_      = 0.0;
  log_delta_range_    = 0.0;

  bin_weight_.assign(1,0.0);
  bin_boundary_.assign(2,min_range_);
  bin_boundary_[1] = max_range_;
}

Distribution::Distribution(const Distribution &dist)
{ // copy constructor
  use_log_spacing_    = dist.use_log_spacing();
  num_bins_           = dist.num_bins();
  min_range_          = dist.min_range();
  max_range_          = dist.max_range();
  total_weight_       = dist.total_weight();
  under_range_weight_ = dist.under_range_weight();
  over_range_weight_  = dist.over_range_weight();
  if(use_log_spacing_) {
    log_delta_range_  = log(max_range_/min_range_); // log(max) - log(min)
    log_min_range_    = log(min_range_);
  } else {
    log_delta_range_ = log_min_range_ = 0.0;
  }

  bin_weight_.assign(num_bins_,0.0);
  bin_boundary_.assign(num_bins_+1,0.0);

  for(int j=0; j<num_bins_; ++j) {
    bin_weight_[j]   = dist.getBinWeight(j);
    bin_boundary_[j] = dist.getBinMin(j);
  }
  bin_boundary_[num_bins_] = dist.getBinMax(num_bins_-1); // max_range_ 
}
Distribution& Distribution::operator=(const Distribution &rdist)
{
  if(this == &rdist) {
    return *this;
  }
  use_log_spacing_    = rdist.use_log_spacing();
  num_bins_           = rdist.num_bins();
  min_range_          = rdist.min_range();
  max_range_          = rdist.max_range();
  total_weight_       = rdist.total_weight();
  under_range_weight_ = rdist.under_range_weight();
  over_range_weight_  = rdist.over_range_weight();
  if(use_log_spacing_) {
    log_delta_range_  = log(max_range_/min_range_); // log(max) - log(min)
    log_min_range_    = log(min_range_);
  } else {
    log_delta_range_ = log_min_range_ = 0.0;
  }
  
  bin_weight_.assign(num_bins_,0.0);
  bin_boundary_.assign(num_bins_+1,0.0);

  for(int j=0; j<num_bins_; ++j) {
    bin_weight_[j]   = rdist.getBinWeight(j);
    bin_boundary_[j] = rdist.getBinMin(j);
  }
  bin_boundary_[num_bins_] = rdist.getBinMax(num_bins_-1); // max_range_ 
  return *this;
}

Distribution& Distribution::operator+=(const Distribution &rdist)
{
  if(isSameBinDistributionAs(rdist)) {
    total_weight_       += rdist.total_weight();
    under_range_weight_ += rdist.under_range_weight();
    over_range_weight_  += rdist.over_range_weight();
    for(int j=0; j<num_bins_; ++j) {
      bin_weight_[j] += rdist.getBinWeight(j);
    }
  } else {
    printf("WARNING: can not use l_dist += r_dist\n");
    printf("                num_bins(): (l) %d   (r) %d\n",
           num_bins_,rdist.num_bins());
    printf("         use_log_spacing(): (l) %s   (r) %s\n",
           (use_log_spacing_ ? "true" : "false"),
           (rdist.use_log_spacing() ? "true" : "false"));
    printf("               min_range(): (l) %.18g   (r) %.18g\n",
           min_range_,rdist.min_range());    
    printf("               max_range(): (l) %.18g   (r) %.18g\n",
           max_range_,rdist.max_range());    
  }
  return *this;
}

Distribution Distribution::operator+(const Distribution &rdist)
{
  Distribution ret_dist; // create a new distribution for return
  if(isSameBinDistributionAs(rdist)) {
    ret_dist=(*this);
    ret_dist+=rdist;
  } else {
    printf("WARNING: can not use l_dist + r_dist\n");
    printf("                num_bins(): (l) %d   (r) %d\n",
           num_bins_,rdist.num_bins());
    printf("         use_log_spacing(): (l) %s   (r) %s\n",
           (use_log_spacing_ ? "true" : "false"),
           (rdist.use_log_spacing() ? "true" : "false"));
    printf("               min_range(): (l) %.18g   (r) %.18g\n",
           min_range_,rdist.min_range());    
    printf("               max_range(): (l) %.18g   (r) %.18g\n",
           max_range_,rdist.max_range());    
  }
  return ret_dist;
}

void Distribution::zeroBinWeights()
{
  total_weight_       = 0.0;
  under_range_weight_ = 0.0;
  over_range_weight_  = 0.0;
  for(int j=0; j<num_bins_; ++j) {
    bin_weight_[j] = 0.0;
  } 
} 

bool Distribution::isSameBinDistributionAs(const Distribution &rdist)
{
  if(num_bins_ != rdist.num_bins()) {
    return false;
  }
  if(min_range_ != rdist.min_range()) {
    return false;
  }
  if(max_range_ != rdist.max_range()) {
    return false;
  }
  if(use_log_spacing_ != rdist.use_log_spacing()) {
    return false;
  }
  return true;
}  

// The linear and logarithmic maps are the direct functions for determining
// the bin index for a particular value.  The map may return an index that is 
// out-of-bounds.  The maps are setup such that if val is in the interval 
// [bin_boundary_[j], bin_boundary_[j+1]) then it is given the index j.  Note 
// that the interval INCLUDES the minimum boundary of the bin and EXCLUDES the
// maximum boundary of the bin.
int Distribution::LinearBinMap(const double val) const
{
  if(isValueInRange(val)) {
    return static_cast<int>(floor((val-min_range_)/(max_range_-min_range_)*
 				  static_cast<double>(num_bins_)));
  } else if(val < min_range_) {
    return -1;
  }
  return num_bins_;
}
int Distribution::LogBinMap(const double val) const
{
  if(isValueInRange(val)) {
    return static_cast<int>(floor((log(val)-log_min_range_)/log_delta_range_*
                                  static_cast<double>(num_bins_)));
  } else if(val < min_range_) {
    return -1;
  }
  return num_bins_;
}

bool Distribution::isBinIndexInRange(const int id) const
{
  if(0 <= id && id < num_bins_) {
    return true;
  }
  return false;
}
bool Distribution::isValueInRange(const double val) const
{
  if(min_range_ <= val && val < max_range_) {
    return true;
  }
  return false;
}
int Distribution::getBinIndex(const double val) const
{
  if(use_log_spacing_) {
    return LogBinMap(val);
  }
  // otherwise use linear map
  return LinearBinMap(val);
}
double Distribution::getBinWeight(const int id) const
{
  if(isBinIndexInRange(id)) {
    return bin_weight_[id];
  }
  else if(id < 0) {
    return under_range_weight_;
  }
  return over_range_weight_;
}
double Distribution::getBinMin(const int id) const
{
  if(isBinIndexInRange(id)) {
    return bin_boundary_[id];
  } 
  printf("WARNING: In Distribution::getBinMin(id),\n");
  printf("         id = %d must be between zero and %d inclusively.\n",
         id,num_bins_-1);
  printf("         Returning NaN.\n");
  return std::numeric_limits<double>::quiet_NaN();
} 

double Distribution::getBinMax(const int id) const
{
  if(isBinIndexInRange(id)) {
    return bin_boundary_[id+1];
  } 
  printf("WARNING: In Distribution::getBinMax(id),\n");
  printf("         id = %d must be between zero and %d inclusively.\n",
         id,num_bins_-1);
  printf("         Returning NaN.\n");
  return std::numeric_limits<double>::quiet_NaN();
}
 
int Distribution::addValue(const double val, const double wt)
{
  int id = getBinIndex(val);
  double previous_total = total_weight_;
  total_weight_ += wt;
  if(total_weight_ == previous_total) {
    printf("WARNING: in Distribution::addValue(value,weight),\n");
    printf("         added weight=%.18g for value=%.18g\n",
           wt,val);
    printf("         does not change the total weight=%.18g\n",
           total_weight_);
  } 
 
  if(isBinIndexInRange(id)) {
    bin_weight_[id]+=wt;
  } else if (id < 0) {
    under_range_weight_+=wt;
  } else {
    over_range_weight_+=wt;
  }
  return id;
}
} // end namespace utilities
} // end namespace zerork
