#include "stdio.h"

#include <string>
#include <iostream>
#include <fstream>

#include "volume_functions.h"

VolumeFromFile::VolumeFromFile(const char filename[], 
                               const InterpolationType interpolation_type,
                               const char *ignore_chars)
{
  int return_flag;
  //printf("Default VolumeFromFile(...) constructor\n"); fflush(stdout);
  return_flag = ReadFile(filename,ignore_chars);
  // TO DO: add handling for a single point to be treated as a constant
  if(return_flag < 2) {
    printf("ERROR: In VolumeFromFile constructor,\n");
    printf("       %d volume time history points read from %s\n",
           return_flag,
           filename);
    fflush(stdout);
    exit(-1);
  }
  num_file_points_ = return_flag;

  SetupInterpolation(interpolation_type);
  volume_multiplier_ = 1.0;
  stroke_multiplier_ = 1.0;
}

VolumeFromFile::VolumeFromFile(const char filename[], 
                               const InterpolationType interpolation_type,
                               const char *ignore_chars,
                               const double volume_multiplier,
                               const double stroke_multiplier)
{
  double min_volume;
  int return_flag;
  //printf("# VolumeFromFile(...) constructor\n"); fflush(stdout);
  //printf("#   with volume_multiplier = %.18g\n",volume_multiplier);
  //printf("#   with stroke_multiplier = %.18g\n",stroke_multiplier);
  fflush(stdout);
  return_flag = ReadFile(filename,ignore_chars);
  // TO DO: add handling for a single point to be treated as a constant
  if(return_flag < 2) {
    printf("ERROR: In VolumeFromFile constructor,\n");
    printf("       %d volume time history points read from %s\n",
           return_flag,
           filename);
    fflush(stdout);
    exit(-1);
  }
  num_file_points_   = return_flag;
  volume_multiplier_ = volume_multiplier;
  stroke_multiplier_ = stroke_multiplier;

  // Rescale the volume with the volume multiplier
  if(volume_multiplier != 1.0) {
    //printf("# Rescaling the volume from the file %s by %.18g\n",
    //       filename,volume_multiplier);
    //fflush(stdout);
    for(int j=0; j<num_file_points_; ++j) {
      volume_[j] *= volume_multiplier;
    }
  }

  // Rescale the volume using the stroke_multiplier:
  // V_new(t) = min(V_file) + (V_file(t) - min(V_file))*stroke_multiplier 
  // The rescaled volume still has the same minimum volume
  if(stroke_multiplier != 1.0) {
    //printf("# Rescaling the stroke from the volume file %s by %.18g\n",
    //       filename,stroke_multiplier);
    //fflush(stdout);
    min_volume = 1.0e300;
    for(int j=0; j<num_file_points_; ++j) {
      if(volume_[j] < min_volume) {
        min_volume = volume_[j];
      }
    }

    for(int j=0; j<num_file_points_; ++j) {
      volume_[j] = min_volume + (volume_[j] - min_volume)*stroke_multiplier;
    }
  }
  // note that the order of applying the volume and stroke_multipliers
  // does not matter
  SetupInterpolation(interpolation_type);
}

int VolumeFromFile::ReadFile(const char filename[],
                             const char *ignore_chars)
{
  size_t ignore_pos;
  std::ifstream input_file(filename, std::ios_base::in);
  std::string line;
  std::string sub_line;
  double time_read, volume_read;
  int num_values_read;
  int num_lines_read = 0;

  t_.clear();
  volume_.clear();

  if(!input_file.is_open()) {
    printf("ERROR: In VolumeFromFile::ReadFile(...),\n");
    printf("       could not open file %s.\n", filename);
    fflush(stdout);
    return -1;
  }

  while(std::getline(input_file, line)) {
    ++num_lines_read;
    // find the position in the string of any comment characters (ignore_chars)
    // beyond which any remaining data is ignored
    ignore_pos = line.find_first_of(ignore_chars);
    sub_line = line.substr(0,ignore_pos);

    if(sub_line.size() > 2) { // minimum length for two values
      num_values_read = sscanf(sub_line.c_str(),
                               "%lf%lf",
                               &time_read,
                               &volume_read);
      if(num_values_read != 2) {
        printf("INFO: Skipping line %d of %s, no data pair read by sscanf:\n",
	       num_lines_read,filename);
        printf("      %s\n",line.c_str());
      } else {
	t_.push_back(time_read);
        volume_.push_back(volume_read);
      }
    } 
  }

  return static_cast<int>(volume_.size());
}

void 
  VolumeFromFile::SetupInterpolation(const InterpolationType interpolation_type)
{
  const bool extrapolate = false; //Make non-const if user settable
  zerork::utilities::InterpolationType interpType = zerork::utilities::LINEAR;
  if(interpolation_type == LINEAR_CLIPPED) {
    //Default
  } else if(interpolation_type == CUBIC_CLIPPED) {
    interpType = zerork::utilities::CUBIC_SPLINE;
  } else {
    printf("ERROR: In VolumeFromFile::SetupInterpolation(),\n");
    printf("       InterpolationType %d not recognized.\n",
           interpolation_type);
    fflush(stdout);
    exit(-1);
  }

  interpolation_table_ = std::make_unique<zerork::utilities::InterpolationTable>(num_file_points_, &t_[0],
                                            &volume_[0], interpType, extrapolate);

  if(interpolation_table_ == NULL) {
    printf("ERROR: In VolumeFromFile::SetupInterpolation(),\n");
    printf("       could not allocate interpolation table.\n");
    exit(-1);
  }  
}

VolumeFromFile::~VolumeFromFile() {}

ReactorError VolumeFromFile::GetVolume(const double reactor_time,
                                       const double state[],
                                       double *volume,
                                       double *dvolume_dt)
{
  *volume = interpolation_table_->Interpolate(reactor_time);
  *dvolume_dt = interpolation_table_->InterpolateFirstDerivative(reactor_time);
  return NONE;
}

