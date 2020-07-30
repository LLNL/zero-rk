#include <stdlib.h>
#include <stdio.h>

#include "integrated_function.h"

IntegratedFunction::IntegratedFunction(double (*function_ptr)(const double t,
                                                              const double y[],
                                                              void *params),
                                       IntegrationType integration_type)
{
  function_ptr_     = function_ptr;
  integration_type_ = integration_type;
  num_points_ = 0;

  if(integration_type_ != TRAPEZOID) {
    printf("WARNING: In IntegratedFunction(...) constructor,\n");
    printf("         integration type %d not recognized.\n",
           integration_type);
    printf("         Setting integration type to TRAPEZOID.\n");
    fflush(stdout);
    integration_type_ = TRAPEZOID;
  }
  // TO DO: set the number of saved points based on the integration type
  num_saved_ = 1;
  // pre-size the previous time and function vectors
  previous_times_.assign(num_saved_,0.0);
  previous_functions_.assign(num_saved_,0.0);

}
double IntegratedFunction::AddPoint(const double t, 
                                    const double y[], 
                                    void *params)
{
  double current_function = function_ptr_(t,y,params);
  if(num_points_ == 0) {
    start_time_       = t;
    current_integral_ = 0.0;
  } else {
    if(t < previous_times_[0]) {
      printf("WARNING: In IntegratedFunction::AddPoint(...),\n");
      printf("         adding point at time t = %.18g,\n", t);
      printf("         that occurs before the last time t = %.18g\n",
             previous_times_[0]);
      fflush(stdout);
    }
    current_integral_ += GetStepArea(t, current_function);
  }
  // update the stored points 
  for(int j = num_saved_-1; j>0; --j) {
    previous_times_[j]     = previous_times_[j-1];
    previous_functions_[j] = previous_functions_[j-1];
  }
  previous_times_[0]     = t;
  previous_functions_[0] = current_function; 
  
  ++num_points_;

  return current_integral_;
}


double IntegratedFunction::GetStepArea(const double t, const double f)
{
  switch(integration_type_) {
    case TRAPEZOID:
      return 0.5*(t-previous_times_[0])*(f+previous_functions_[0]);

    default:
      printf("WARNING: In IntegratedFunction(...) constructor,\n");
      printf("         integration type %d not recognized.\n",
             integration_type_);
      printf("         Setting integration type to TRAPEZOID.\n");
      fflush(stdout);
      integration_type_ = TRAPEZOID;
      return GetStepArea(t,f);
  }
}


double IntegratedFunction::GetCurrentIntegral() const
{
  if(num_points_ < 1) {
    printf("WARNING: In IntegratedFunction::GetCurrentIntegral(),\n");
    printf("         the number of points added is %d.\n",num_points_);
    printf("         Returning zero.\n");
    fflush(stdout);
    return 0.0;
  }
  return current_integral_;
}
double IntegratedFunction::GetStartTime() const
{
  if(num_points_ < 1) {
    printf("WARNING: In IntegratedFunction::GetStartTime(),\n");
    printf("         the number of points added is %d.\n",num_points_);
    printf("         Returning zero.\n");
    fflush(stdout);
    return 0.0;
  }
  return start_time_;
}  
double IntegratedFunction::GetLastTime() const
{
  if(num_points_ < 1) {
    printf("WARNING: In IntegratedFunction::GetLastTime(),\n");
    printf("         the number of points added is %d.\n",num_points_);
    printf("         Returning zero.\n");
    fflush(stdout);
    return 0.0;
  }
  return previous_times_[0];
}
double IntegratedFunction::GetLastFunction() const
{
  if(num_points_ < 1) {
    printf("WARNING: In IntegratedFunction::GetLastFunction(),\n");
    printf("         the number of points added is %d.\n",num_points_);
    printf("         Returning zero.\n");
    fflush(stdout);
    return 0.0;
  }
  return previous_functions_[0];
}
