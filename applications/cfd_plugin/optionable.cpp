
#include "optionable.h"

int Optionable::GetIntOption(const std::string option_name, int* option_value) {
  std::map<std::string,int>::iterator it;
  it = int_options_.find(option_name);
  if( it != int_options_.end() ) {
    *option_value = it->second;
    return 0;
  } else {
    return 1;
  }
}

int Optionable::GetDoubleOption(const std::string option_name, double* option_value) {
  std::map<std::string,double>::iterator it;
  it = double_options_.find(option_name);
  if( it != double_options_.end() ) {
    *option_value = it->second;
    return 0;
  } else {
    return 1;
  }
}

void Optionable::GetIntOptions(IntOptions* options) {
  *options = int_options_;
}

void Optionable::GetDoubleOptions(DoubleOptions* options) {
  *options = double_options_;
}

int Optionable::SetIntOption(const std::string option_name, int option_value) {
  int_options_[option_name] = option_value;
  return 0;
}

int Optionable::SetDoubleOption(const std::string option_name, double option_value) {
  double_options_[option_name] = option_value;
  return 0;
}

void Optionable::SetIntOptions(const IntOptions& options) {
  int_options_ = options;
}

void Optionable::SetDoubleOptions(const DoubleOptions& options) {
  double_options_ = options;
}


