#ifndef OPTIONABLE_H
#define OPTIONABLE_H

#include "optionable_base.h"

class Optionable : public OptionableBase
{
 public:
  Optionable() {};
  ~Optionable() {};

  int GetIntOption(const std::string option_name, int* option_value);
  int GetDoubleOption(const std::string option_name, double* option_value);

  void GetIntOptions(IntOptions* options);
  void GetDoubleOptions(DoubleOptions* options);

  int SetIntOption(const std::string option_name, int option_value);
  int SetDoubleOption(const std::string option_name, double option_value);

  void SetIntOptions(const IntOptions& options);
  void SetDoubleOptions(const DoubleOptions& options);
};


#endif
