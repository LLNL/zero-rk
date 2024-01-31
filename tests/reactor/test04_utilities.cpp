#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "test04_utilities.h"

int GetForwardFirstOrder(JacobianRhsFunction function,
                         const double relative_delta,
                         const double absolute_delta,
                         const size_t num_equations,
                         const double t,
                         const double y[],
                         void *params,
                         double jacobian[],
                         double truncation_error[],
                         double roundoff_error[]);

int GetCentralFourthOrder(JacobianRhsFunction function,
                          const double relative_delta,
                          const double absolute_delta,
                          const size_t num_equations,
                          const double t,
                          const double y[],
                          void *params,
                          double jacobian[],
                          double truncation_error[],
                          double roundoff_error[]);

void CentralDerivativeFourthOrder(const double f_stencil[],
                                  const double y,
                                  const double dy,
                                  double *result,
                                  double *absolute_truncation,
                                  double *absolute_roundoff);
                                   

int GetJacobian(const JacobianMethod method,
                JacobianRhsFunction function,
                const double relative_delta,
                const double absolute_delta,
                const size_t num_equations,
                const double t,
                const double y[],
                void *params,
                double jacobian[],
                double truncation_error[],
                double roundoff_error[])
{
  if(method == FORWARD_FIRST_ORDER) {
    
    return GetForwardFirstOrder(function,
                                relative_delta,
                                absolute_delta,
                                num_equations,
                                t,
                                y,
                                params,
                                jacobian,
                                truncation_error,
                                roundoff_error);
  }
  // if CENTRAL_FOURTH_ORDER
  if(method == CENTRAL_FOURTH_ORDER) {
    
    return GetCentralFourthOrder(function,
                                 relative_delta,
                                 absolute_delta,
                                 num_equations,
                                 t,
                                 y,
                                 params,
                                 jacobian,
                                 truncation_error,
                                 roundoff_error);
  }

  printf("WARNING: In GetJacobian(...),\n");
  printf("         method type %d (enum) not recognized.\n",
         method);
    
  return -1;
}

int GetForwardFirstOrder(JacobianRhsFunction function,
                         const double relative_delta,
                         const double absolute_delta,
                         const size_t num_equations,
                         const double t,
                         const double y[],
                         void *params,
                         double jacobian[],
                         double truncation_error[],
                         double roundoff_error[])
{
  const double machine_epsilon = 2.2204460492503131e-16;
  double delta_y, inv_delta_y;
  double *perturb_y, *perturb_f, *f;
  int flag;

  f         = new double[num_equations];

  flag = function(num_equations, t, y, f, params);
  if(flag != 0) {
    delete [] f;
    return flag;
  }

  perturb_f = new double[num_equations];
  perturb_y = new double[num_equations];

  for(size_t j=0; j<num_equations; ++j) {

    delta_y = fabs(y[j]*relative_delta) + fabs(absolute_delta);

    for(size_t k=0; k<num_equations; ++k) {
      perturb_y[k] = y[k];
    } 
    perturb_y[j] += delta_y;
    // recompute delta from the finite precision representation of y+dy and y
    delta_y = perturb_y[j] - y[j];
    inv_delta_y = 1.0/delta_y;

    flag = function(num_equations, t, perturb_y, perturb_f, params);
    if(flag != 0) {
      delete [] f;
      delete [] perturb_f;
      delete [] perturb_y;
      return flag;
    }

    for(size_t k=0; k<num_equations; ++k) {

      jacobian[k+j*num_equations] = (perturb_f[k] - f[k])*inv_delta_y;
      truncation_error[k+j*num_equations] = 1.0e300; // no formula 
      roundoff_error[k+j*num_equations] = 
        fabs(f[k]*machine_epsilon*inv_delta_y);
    }
  }

  delete [] f;
  delete [] perturb_f;
  delete [] perturb_y;
  return flag;
}

int GetCentralFourthOrder(JacobianRhsFunction function,
                          const double relative_delta,
                          const double absolute_delta,
                          const size_t num_equations,
                          const double t,
                          const double y[],
                          void *params,
                          double jacobian[],
                          double truncation_error[],
                          double roundoff_error[])
{
  double initial_result, optimum_result;
  double initial_delta_y, optimum_delta_y;
  double initial_truncation_error, optimum_truncation_error;
  double initial_roundoff_error, optimum_roundoff_error;
  double initial_total_error, optimum_total_error;
  double *perturb_y, *perturb_f;
  double f_stencil[4];
  int flag = 0;
  
  perturb_f = new double[num_equations];
  perturb_y = new double[num_equations];

  for(size_t j=0; j<num_equations; ++j) { // column j

    for(size_t k=0; k<num_equations; ++k) { // row k
      
      // copy y to the perturbation vector
      for(size_t m=0; m<num_equations; ++m) {
        perturb_y[m] = y[m];
      }
      // compute initial perturbation for y[j]
      initial_delta_y = fabs(relative_delta*perturb_y[j]) + 
        fabs(absolute_delta);
      perturb_y[j] += initial_delta_y;
      initial_delta_y = perturb_y[j] - y[j];

      // f_stencil[3] = f(y + dy)
      flag=function(num_equations,
                    t,
                    perturb_y,
                    perturb_f,
                    params);
      if(flag != 0) {
        delete [] perturb_f;
        delete [] perturb_y;
        return flag;
      }
      f_stencil[3] = perturb_f[k];

      // f_stencil[2] = f(y + dy/2)
      perturb_y[j] = y[j] + initial_delta_y/2.0;
      flag=function(num_equations,
                    t,
                    perturb_y,
                    perturb_f,
                    params);
      if(flag != 0) {
        delete [] perturb_f;
        delete [] perturb_y;
        return flag;
      }
      f_stencil[2] = perturb_f[k];

      // f_stencil[1] = f(y - dy/2)
      perturb_y[j] = y[j] - initial_delta_y/2.0;
      flag=function(num_equations,
                    t,
                    perturb_y,
                    perturb_f,
                    params);
      if(flag != 0) {
        delete [] perturb_f;
        delete [] perturb_y;
        return flag;
      }
      f_stencil[1] = perturb_f[k];


      // f_stencil[0] = f(y - dy)
      perturb_y[j] = y[j] - initial_delta_y;
      flag=function(num_equations,
                    t,
                    perturb_y,
                    perturb_f,
                    params);
      if(flag != 0) {
        delete [] perturb_f;
        delete [] perturb_y;
        return flag;
      }
      f_stencil[0] = perturb_f[k];

      CentralDerivativeFourthOrder(f_stencil,
                                   y[j],
                                   initial_delta_y,
                                   &initial_result,
                                   &initial_truncation_error,
                                   &initial_roundoff_error);

      initial_total_error = initial_truncation_error+initial_roundoff_error;

      if((initial_roundoff_error < initial_truncation_error) &&
         (initial_roundoff_error > 0.0 && initial_truncation_error > 0.0)) {

        // compute the optimum delta_y to to minimize the total error
        optimum_delta_y = 
          initial_delta_y*pow(initial_roundoff_error/
                              (2.0*initial_truncation_error),
                              1.0/3.0);
        // f_stencil[3] = f(y + dy)
        perturb_y[j] = y[j] + optimum_delta_y;
        flag=function(num_equations,
                      t,
                      perturb_y,
                      perturb_f,
                      params);
        if(flag != 0) {
          delete [] perturb_f;
          delete [] perturb_y;
          return flag;
        }
        f_stencil[3] = perturb_f[k];

        // f_stencil[2] = f(y + dy/2)
        perturb_y[j] = y[j] + optimum_delta_y/2.0;
        flag=function(num_equations,
                      t,
                      perturb_y,
                      perturb_f,
                      params);
        if(flag != 0) {
          delete [] perturb_f;
          delete [] perturb_y;
          return flag;
        }
        f_stencil[2] = perturb_f[k];

        // f_stencil[1] = f(y - dy/2)
        perturb_y[j] = y[j] - optimum_delta_y/2.0;
        flag=function(num_equations,
                      t,
                      perturb_y,
                      perturb_f,
                      params);
        if(flag != 0) {
          delete [] perturb_f;
          delete [] perturb_y;
          return flag;
        }
        f_stencil[1] = perturb_f[k];

        // f_stencil[0] = f(y - dy)
        perturb_y[j] = y[j] - optimum_delta_y;
        flag=function(num_equations,
                      t,
                      perturb_y,
                      perturb_f,
                      params);
        if(flag != 0) {
          delete [] perturb_f;
          delete [] perturb_y;
          return flag;
        }
        f_stencil[0] = perturb_f[k];

        CentralDerivativeFourthOrder(f_stencil,
                                     y[j],
                                     optimum_delta_y,
                                     &optimum_result,
                                     &optimum_truncation_error,
                                     &optimum_roundoff_error);

        optimum_total_error = optimum_truncation_error+optimum_roundoff_error;

        // overwrite the initial result and error if optimum is considered
        // a valid improvement
        if(optimum_total_error < initial_total_error &&
           fabs(optimum_result - initial_result) < 4.0*initial_total_error) {
          
          initial_result           = optimum_result;
          initial_truncation_error = optimum_truncation_error;
          initial_roundoff_error   = optimum_roundoff_error; 
        }

      } // if optimum_delta_y can be calculated 
      
      jacobian[k+j*num_equations]         = initial_result;
      truncation_error[k+j*num_equations] = initial_truncation_error;
      roundoff_error[k+j*num_equations]   = initial_roundoff_error ;

    } // for loop over row k
  } // for loop over column j

  //flag = function(num_equations, t, y, f, params);
  if(flag != 0) {

    delete [] perturb_f;
    delete [] perturb_y;
    return flag;
  }

  delete [] perturb_f;
  delete [] perturb_y;

  return flag;
}

// f_stencil[0] = f(y - dy)
// f_stencil[1] = f(y - dy/2)
// f_stencil[2] = f(y + dy/2)
// f_stencil[3] = f(y + dy)
void CentralDerivativeFourthOrder(const double f_stencil[],
                                  const double y,
                                  const double dy,
                                  double *result,
                                  double *absolute_truncation,
                                  double *absolute_roundoff)
{
  const double machine_epsilon = 2.2204460492503131e-16;
  double result_order2 = 0.5*(f_stencil[3]-f_stencil[0]); // r3 in gsl
  double result_order4 = (4.0/3.0)*(f_stencil[2]-f_stencil[1]) - 
    (1.0/3.0)*result_order2; // r5 in gsl
  double error_order2 = 
    (fabs(f_stencil[3])+fabs(f_stencil[0]))*machine_epsilon; // e3 in gsl
  double error_order4 = 
    2.0*(fabs(f_stencil[2])+fabs(f_stencil[1]))*machine_epsilon + error_order2;
  double finite_precision = 
    ((fabs(result_order2/dy) > fabs(result_order4/dy)) ?
     fabs(result_order2/dy) :
     fabs(result_order4/dy));

  finite_precision *= (fabs(y)/dy)*machine_epsilon;
  
  *result              = result_order4/dy;
  *absolute_truncation = fabs((result_order4-result_order2)/dy);
  *absolute_roundoff   = fabs(error_order4/dy) + finite_precision;
}
