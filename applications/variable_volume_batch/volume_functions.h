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

  VolumeFromFile(const char filename[], 
                 const InterpolationType interpolation_type,
                 const char *ignore_chars,
                 const double volume_multiplier,
                 const double stroke_multiplier,
                 const bool is_csv,
                 const int num_lines_skipped);


  ~VolumeFromFile();
  ReactorError GetVolume(const double reactor_time,
                         const double state[],
                         double *volume,
                         double *dvolume_dt);
  double GetMinTime() const {return t_[0];}
  double GetMaxTime() const {return t_[num_file_points_-1];}
  int GetNumFilePoints() const {return num_file_points_;}

  double GetMaxVolume() const {return max_volume_;}
  double GetMinVolume() const {return min_volume_;}
  double GetTimeAtMaxVolume() const {return time_at_max_volume_;}
  double GetTimeAtMinVolume() const {return time_at_min_volume_;}


 private:
  int ReadFile(const char filename[], const char *ignore_chars);
  int ReadCSVFile(const char filename[], 
                  const char *ignore_chars,
                  const int num_lines_skipped);
  void SetFileExtrema();
  void SetupInterpolation(const InterpolationType interpolation_type);
  int num_file_points_;
  double volume_multiplier_;
  double stroke_multiplier_;

  // file extrema
  double max_volume_;
  double min_volume_;
  double time_at_min_volume_;
  double time_at_max_volume_;

  std::vector<double> t_;
  std::vector<double> volume_;

  std::unique_ptr<zerork::utilities::InterpolationTable> interpolation_table_;
};

#endif
