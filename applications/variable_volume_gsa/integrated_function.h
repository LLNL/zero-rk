#ifndef INTEGRATED_FUNCTION_H_
#define INTEGRATED_FUNCTION_H_

#include <vector>

enum IntegrationType {TRAPEZOID};

class IntegratedFunction
{
 public:
  IntegratedFunction(double (*function_ptr)(const double t,
                                            const double y[],
                                            void *params),
                     IntegrationType integration_type);
  double AddPoint(const double t, const double y[], void *params);

  int GetNumPoints() const {return num_points_;}
  double GetCurrentIntegral() const;
  double GetStartTime() const;  
  double GetLastTime() const;
  double GetLastFunction() const;

 private:
  double GetStepArea(const double t, const double f);
  
  IntegrationType integration_type_;
  int num_points_;
  int num_saved_;
  double (*function_ptr_)(const double t,
                          const double y[],
                          void *params);
  std::vector<double> previous_times_;
  std::vector<double> previous_functions_;
  double current_integral_;
  double start_time_;
};


#endif
