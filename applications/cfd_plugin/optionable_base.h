#ifndef OPTIONABLE_BASE_H
#define OPTIONABLE_BASE_H

#include <string>
#include <map>

class OptionableBase
{
 public:
  OptionableBase() {};
  virtual ~OptionableBase() {};

  typedef std::map<std::string,int> IntOptions;
  typedef std::map<std::string,double> DoubleOptions;
  typedef std::map<std::string,std::string> StringOptions;

  virtual int GetIntOption(const std::string option_name, int* option_value) = 0;
  virtual int GetDoubleOption(const std::string option_name, double* option_value) = 0;
  virtual int GetStringOption(const std::string option_name, std::string* option_value) = 0;

  virtual void GetIntOptions(IntOptions* options) = 0;
  virtual void GetDoubleOptions(DoubleOptions* options) = 0;
  virtual void GetStringOptions(StringOptions* options) = 0;

  virtual int SetIntOption(const std::string option_name, int option_value) = 0;
  virtual int SetDoubleOption(const std::string option_name, double option_value) = 0;
  virtual int SetStringOption(const std::string option_name, std::string option_value) = 0;

  virtual void SetIntOptions(const IntOptions& options) = 0;
  virtual void SetDoubleOptions(const DoubleOptions& options) = 0;
  virtual void SetStringOptions(const StringOptions& options) = 0;

 protected:
  IntOptions int_options_;
  DoubleOptions double_options_;
  StringOptions string_options_;
};


#endif
