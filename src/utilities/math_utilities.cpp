#include <stdlib.h>
#include <stdio.h>
#ifdef _WIN32
#define _USE_MATH_DEFINES //for M_PI
#endif
#include <cmath>

#include <map>

#include "math_utilities.h"

namespace zerork
{
namespace utilities
{


#ifdef _WIN32
//The following drand48 implementation is courtesy:
// https://gist.github.com/mortennobel/8665258
// (accessed 20220407)
// Originally from FreeBSD:
/*
 * Copyright (c) 1993 Martin Birgmeier
 * All rights reserved.
 *
 * You may redistribute unmodified or modified versions of this source
 * code provided that the above copyright notice and this and the
 * following conditions are retained.
 *
 * This software is provided ``as is'', and comes with no warranties
 * of any kind. I shall in no event be liable for anything that happens
 * to anyone/anything when using this software.
 */
#define RAND48_SEED_0 (0x330e)
#define RAND48_SEED_1 (0xabcd)
#define RAND48_SEED_2 (0x1234)
#define RAND48_MULT_0 (0xe66d)
#define RAND48_MULT_1 (0xdeec)
#define RAND48_MULT_2 (0x0005)
#define RAND48_ADD (0x000b)

unsigned short _rand48_seed[3] = {
        RAND48_SEED_0,
         RAND48_SEED_1,
         RAND48_SEED_2
};
unsigned short _rand48_mult[3] = {
         RAND48_MULT_0,
         RAND48_MULT_1,
         RAND48_MULT_2
 };
unsigned short _rand48_add = RAND48_ADD;

void
 _dorand48(unsigned short xseed[3])
 {
	         unsigned long accu;
	         unsigned short temp[2];
	
	         accu = (unsigned long)_rand48_mult[0] * (unsigned long)xseed[0] +
	          (unsigned long)_rand48_add;
	         temp[0] = (unsigned short)accu;        /* lower 16 bits */
	         accu >>= sizeof(unsigned short)* 8;
	         accu += (unsigned long)_rand48_mult[0] * (unsigned long)xseed[1] +
	          (unsigned long)_rand48_mult[1] * (unsigned long)xseed[0];
	         temp[1] = (unsigned short)accu;        /* middle 16 bits */
	         accu >>= sizeof(unsigned short)* 8;
	         accu += _rand48_mult[0] * xseed[2] + _rand48_mult[1] * xseed[1] + _rand48_mult[2] * xseed[0];
	         xseed[0] = temp[0];
	         xseed[1] = temp[1];
	         xseed[2] = (unsigned short)accu;
}

double erand48(unsigned short xseed[3])
{
         _dorand48(xseed);
         return ldexp((double) xseed[0], -48) +
                ldexp((double) xseed[1], -32) +
                ldexp((double) xseed[2], -16);
}

double drand48(){
	return erand48(_rand48_seed);
}

void srand48(long seed){
	_rand48_seed[0] = RAND48_SEED_0;
	_rand48_seed[1] = (unsigned short)seed;
	_rand48_seed[2] = (unsigned short)(seed >> 16);
	_rand48_mult[0] = RAND48_MULT_0;
	_rand48_mult[1] = RAND48_MULT_1;
	_rand48_mult[2] = RAND48_MULT_2;
	_rand48_add = RAND48_ADD;
}
#endif

double random01() {
   return drand48();
}

void random01seed(int seed) {
   srand48(seed);
}

void SortVectors(std::vector<double> *keys,
                 std::vector<double> *values)
{
  const size_t num_elements = ((keys->size() < values->size()) ?
                               keys->size() : values->size());
  std::multimap<double, double> sort_map;
  // Need typename because std::multimap<KEY, VALUE> is dependent scope
  typename std::multimap<double, double>::iterator iter;

  // The iterator of the last insertion can be used as a hint to improve
  // efficiency when the list is already sorted.
  // iter = sort_map.begin(); // set first hint
  for(size_t j=0; j<num_elements; ++j) {
    sort_map.insert(std::pair<double, double>(keys->at(j), values->at(j)));
    // use hint of the last position
    //iter = sort_map.insert(iter, keys[j], values[j]);
  }
  // copy the map to key and value arrays
  size_t k=0;
  for(iter=sort_map.begin(); iter != sort_map.end(); ++iter) {
    keys->at(k)   = iter->first;
    values->at(k) = iter->second;
    ++k;
  }

}

// template<typename KEY, typename VALUE>
// void SortVectors(std::vector<KEY> &keys,
//                  std::vector<VALUE> &values)
// {
//   const size_t num_elements = ((keys.size() < values.size()) ?
//                                keys.size() : values.size());
//   std::multimap<KEY, VALUE> sort_map;
//   // Need typename because std::multimap<KEY, VALUE> is dependent scope
//   typename std::multimap<KEY, VALUE>::iterator iter;

//   // The iterator of the last insertion can be used as a hint to improve
//   // efficiency when the list is already sorted.
//   // iter = sort_map.begin(); // set first hint
//   for(size_t j=0; j<num_elements; ++j) {
//     sort_map.insert(std::pair<KEY, VALUE>(keys[j], values[j]));
//     // use hint of the last position
//     //iter = sort_map.insert(iter, keys[j], values[j]);
//   }
//   // copy the map to key and value arrays
//   size_t k=0;
//   for(iter=sort_map.begin(); iter != sort_map.end(); ++iter) {
//     keys[k]   = iter->first;
//     values[k] = iter->second;
//     ++k;
//   }
// }

InterpolationTable::InterpolationTable(const int num_points,
                          const double *x,
			  const double *f,
                          const InterpolationType interpolation_type,
                          const bool use_extrapolation)
{
  if(num_points < 2) {
    // create a default table
    num_points_ = 2;
    x_.push_back(-1.0e300);
    x_.push_back(1.0e300);
    f_.assign(num_points,0.0);

    printf("# WARNING: In InterpolationTable constructor,\n");
    printf("#          number of table points %d can not be used.\n",
           num_points);
    printf("#          Creating default two point table:\n");
    printf("#              f(x[0]) = %24.16e, x[0] = %24.16e\n",f_[0],x_[0]);
    printf("#              f(x[1]) = %24.16e, x[1] = %24.16e\n",f_[1],x_[1]);
  } else {

    std::vector<double> x_copy, f_copy;
    x_copy.assign(num_points, 0.0);
    f_copy.assign(num_points, 0.0);
    for(int j=0; j<num_points; ++j) {
      x_copy[j] = x[j];
      f_copy[j] = f[j];
    }
    SortVectors(&x_copy, &f_copy);

    x_.push_back(x_copy[0]);
    f_.push_back(f_copy[0]);

    for(int j=1; j<num_points; ++j) {
      // only record points that are strictly increasing
      if(x_copy[j] > x_copy[j-1]) {
        x_.push_back(x_copy[j]);
        f_.push_back(f_copy[j]);
      } else {
        printf("# WARNING: In InterpolationTable constructor,\n");
        printf("#          skipping non-increasing table point:\n");
        printf("#              f(x[%d]) = %24.16e, x[%d] = %24.16e\n",
               j, f_copy[j], j, x_copy[j]);
      }
    }
    num_points_ = (int)x_.size();
  }


  last_id_ = 0;

  interpolation_type_ = interpolation_type;
  use_extrapolation_ = use_extrapolation;

  SetWorkspace();
}

void InterpolationTable::SetWorkspace()
{
  const int num_intervals = num_points_-1;

  // TODO: add other interpolation types
  if(interpolation_type_ == CUBIC_SPLINE) {
    //Implementing algorithm for natural cubic spline interpolation:
    //  Numerical Analysis 9th ed., R.L. Burden & J.D. Faires 
    //  pp. 149-150
    int num_coeffs = 3*num_points_;
    workspace_.assign(num_coeffs, 0.0); //b, c, d
    std::vector<double> h(num_intervals);
    std::vector<double> alpha(num_intervals);
    std::vector<double> mu(num_intervals);
    std::vector<double> l(num_points_);
    std::vector<double> z(num_points_);

    for(int j=0; j<num_intervals; ++j) {
      h[j] = (x_[j+1]-x_[j]);
    }
    for(int j=1; j<num_intervals; ++j) {
      alpha[j] = 3 * (f_[j+1]-f_[j]) / h[j] - 3 * (f_[j]-f_[j-1]) / h[j-1];
    }

    l[0] =  1;
    mu[0] = 0;
    z[0] =  0;

    for(int j=1; j<num_intervals; ++j) {
      l[j] = 2*(x_[j+1]-x_[j-1]) - h[j-1]*mu[j-1];
      mu[j] = h[j]/l[j];
      z[j] = (alpha[j] - h[j-1]*z[j-1])/l[j];
    }

    l[num_points_-1] = 1;
    z[num_points_-1] = 0;
    workspace_[num_intervals*3+1] = 0; //c_n

    for(int j=num_intervals-1; j>= 0; --j) {
      double c_jp1 = workspace_[(j+1)*3+1];
      double c_j   = z[j] - mu[j]*c_jp1;
      workspace_[j*3+0] = (f_[j+1] - f_[j])/h[j]
                          - h[j]*(c_jp1 + 2*c_j)/3; //b_j
      workspace_[j*3+1] = c_j; //c_j
      workspace_[j*3+2] = (c_jp1-c_j)/(3*h[j]); //d_j
    }
  } else {
    // interpolation_type_ == LINEAR is default
    workspace_.assign(num_intervals, 0.0);
    for(int j=0; j<num_intervals; ++j) {
      workspace_[j] = (f_[j+1]-f_[j])/(x_[j+1]-x_[j]);
    }
  }
}

int InterpolationTable::GetIntervalId(const double x) const
{
  const int num_intervals = num_points_-1;

  if(x < x_[0]) {
    return -1;
  } else if(x > x_[num_points_-1]) {
    return num_intervals;
  }

  //// check if last_id_ interval is still valid
  //if(0 <= last_id_ && last_id_ <= num_points_-2) {
  //
  //  if(x_[last_id_] <= x && x < x_[last_id_+1]) {
  //    return last_id_;
  //  }
  //}
  // TODO: replace with a nice binary search
  for(int j=last_id_; j<num_intervals; ++j) {
    if(x_[j] <= x && x <= x_[j+1]) {
      last_id_ = j;
      return j;
    }
  }
  for(int j=0; j<last_id_; ++j) {
    if(x_[j] <= x && x <= x_[j+1]) {
      last_id_ = j;
      return j;
    }
  }
  printf("# ERROR: could not find interval for point x=%24.18e\n",
         x);
  exit(-1);
  return 0;
}

double InterpolationTable::Interpolate(const double x) const
{
  int id;
  double f_interp;
  // deterimine the interval id of the value x

  if(x_[0] < x && x < x_[num_points_-1]) {
    id = GetIntervalId(x);
  } else if(x <= x_[0]) {

    if(use_extrapolation_) {
      id = 0; // compute f using the first interval
    } else {
      return f_[0];
    }

  } else {

    if(use_extrapolation_) {
      id = num_points_-2; // compute f using the last interval
    } else {
      return f_[num_points_-1];
    }

  }

  // TODO: add other interpolation types
  if(interpolation_type_ == CUBIC_SPLINE) {
    const double dx = x-x_[id];
    const double bi = workspace_[3*id];
    const double ci = workspace_[3*id+1];
    const double di = workspace_[3*id+2];
    f_interp = f_[id] + dx*(bi + dx*(ci + dx*di));
  } else {
  // interpolation_type_ == LINEAR is default
    f_interp = f_[id] + workspace_[id]*(x - x_[id]);
  }

  return f_interp;
}


double InterpolationTable::InterpolateFirstDerivative(const double x) const
{
  int id;
  double df_interp;
  // deterimine the interval id of the value x

  if(x_[0] < x && x < x_[num_points_-1]) {
    id = GetIntervalId(x);
  } else if(x <= x_[0]) {
    if(use_extrapolation_) {
      id = 0; // compute f using the first interval
    } else {
      return 0.0;
    }
  } else {
    if(use_extrapolation_) {
      id = num_points_-2; // compute f using the last interval
    } else {
      return 0.0;
    }
  }

  if(interpolation_type_ == CUBIC_SPLINE) {
    const double dx = x-x_[id];
    const double bi = workspace_[3*id];
    const double ci = workspace_[3*id+1];
    const double di = workspace_[3*id+2];
    df_interp = bi + dx*(2*ci + 3*di*dx);
  } else {
  // interpolation_type_ == LINEAR is default
    df_interp = workspace_[id];
  }

  return df_interp;
}

// Approximate inverse error function implemented from:
// "A handy approximation for the error function and its inverse"
// by Sergei Winitzki
// Note we are using the erf_inv form and transforming q to p
// Worst case error should be less than 1% over whole range
template<typename T>
T erfc_inv(T q) {
    if(q <= 0) return INFINITY;
    if(q >= 2) return -INFINITY;
    T p = 1 - q;
    T sgn = p < 0 ? -1 : 1;

    T logpt = log((1 - p)*(1 + p));
    T pt1 = 2/(M_PI*0.147) + 0.5*logpt;
    T pt2 = 1/(0.147) * logpt;

    return(sgn*sqrt(-pt1 + sqrt(pt1*pt1 - pt2)));
}

template float erfc_inv<float>(float q);
template double erfc_inv<double>(double q);
template long double erfc_inv<long double>(long double q);

} // namespace utilities
} // namespace zerork
