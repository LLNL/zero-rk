#ifndef MATH_UTILITIES_H_
#define MATH_UTILITIES_H_

#include <vector>

namespace zerork 
{
namespace utilities
{

void SortVectors(std::vector<double> *keys, 
                 std::vector<double> *values);

// TODO: Implement templated version
//
// SortVectors sorts the first num_elements of the value array using the
// associated key array.  The key and value arrays can have different types
// specified by the template.  The KEY type needs to support a comparison
// operation to allow comparison for ordering purposes (e.g., char, int, float,
// and double).
//
// Inputs:
//     num_elements:  number of elements to sort
//
// Inputs & Outputs:
//    keys:           array of keys (of type KEY) used to sort the 
//                    associated array in ascending order. On return, the
//                    keys arrays will be in ascending order.
//    values:         array of values (of type VALUE) that is sorted by
//                    the associsiated array of keys in ascending order. This
//                    means that on return values[0] will have the smallest
//                    key in the keys array, values[1] will have the next
//                    largest and then so on.
//
// Notes:
//
// (1) There are currently no checks that the type KEY is sortable.
// (2) The values and keys are inserted in the std::multimap container so
//     the order of values with the same key depends on the operation of the
//     multimap insert function
//template<typename KEY, typename VALUE> 
//  void SortVectors(std::vector<KEY> &keys, 
//                   std::vector<VALUE> &values);



enum InterpolationType {LINEAR, CUBIC_SPLINE};

class InterpolationTable
{
 public:
  InterpolationTable(const int num_points,
                     const double *x,
                     const double *f,
                     const InterpolationType interpolation_type,
                     const bool use_extrapolation);
  ~InterpolationTable() {};
  double Interpolate(const double x) const;
  double InterpolateFirstDerivative(const double x) const;

 private:
  void SetWorkspace();
  int GetIntervalId(const double x) const;

  mutable int last_id_;
  
  bool use_extrapolation_;
  InterpolationType interpolation_type_;
  int num_points_;
  std::vector<double> x_;
  std::vector<double> f_;
  std::vector<double> workspace_;
};

template<typename T>
T erfc_inv(T q);

double random01();

void random01seed(int seed);

} // namespace utilities
} // namespace advcomb

#endif
