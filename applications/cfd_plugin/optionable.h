#ifndef OPTIONABLE_H
#define OPTIONABLE_H

#include "optionable_base.h"

class Optionable : public OptionableBase
{
 public:
  Optionable() {};
  ~Optionable() {};

  optionable_status_t GetIntOption(const std::string option_name, int* option_value);
  optionable_status_t GetDoubleOption(const std::string option_name, double* option_value);
  optionable_status_t GetStringOption(const std::string option_name, std::string* option_value);

  void GetIntOptions(IntOptions* options);
  void GetDoubleOptions(DoubleOptions* options);
  void GetStringOptions(StringOptions* options);

  optionable_status_t SetIntOption(const std::string option_name, int option_value);
  optionable_status_t SetDoubleOption(const std::string option_name, double option_value);
  optionable_status_t SetStringOption(const std::string option_name, std::string option_value);

  void SetIntOptions(const IntOptions& options);
  void SetDoubleOptions(const DoubleOptions& options);
  void SetStringOptions(const StringOptions& options);
};


#endif
