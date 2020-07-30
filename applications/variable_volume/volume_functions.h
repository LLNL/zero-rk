#ifndef VOLUME_FUNCTIONS_H_
#define VOLUME_FUNCTIONS_H_

#include <memory>
#include <vector>

#include "math_utilities.h"
#include "reactor/reactor_constants.h"

enum InterpolationType {LINEAR_CLIPPED, CUBIC_CLIPPED};

class VolumeFromFile : public Volume
{
 public:
  VolumeFromFile(const char filename[], 
                 const InterpolationType interpolation_type,
                 const char *ignore_chars);
  VolumeFromFile(const char filename[], 
                 const InterpolationType interpolation_type,
                 const char *ignore_chars,
                 const double volume_multiplier,
                 const double stroke_multiplier);

  ~VolumeFromFile();
  ReactorError GetVolume(const double reactor_time,
                         const double state[],
                         double *volume,
                         double *dvolume_dt);
  double GetMinTime() const {return t_[0];}
  double GetMaxTime() const {return t_[num_file_points_-1];}
  int GetNumFilePoints() const {return num_file_points_;}

 private:
  int ReadFile(const char filename[], const char *ignore_chars);
  void SetupInterpolation(const InterpolationType interpolation_type);
  int num_file_points_;
  double volume_multiplier_;
  double stroke_multiplier_;

  std::vector<double> t_;
  std::vector<double> volume_;

  std::unique_ptr<zerork::utilities::InterpolationTable> interpolation_table_;
};

#endif
