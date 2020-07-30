#ifndef TEST_UTILITIES_H_
#define TEST_UTILITIES_H_

enum JacobianMethod {FORWARD_FIRST_ORDER, CENTRAL_FOURTH_ORDER};

typedef int (*JacobianRhsFunction)(const size_t num_equations,
                                   const double t, 
                                   const double y[],
                                   double f[],
                                   void *params) ;


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
                double roundoff_error[]);



#endif
