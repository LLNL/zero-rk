#ifndef ZERORK_EXTERNAL_FUNCS_H
#define ZERORK_EXTERNAL_FUNCS_H

namespace zerork {

typedef bool (*external_func_check_t)(const int,const int);
typedef void (*external_func_rates_t)(const double [], double [], double [], double []);
typedef void (*external_func_arrh_t)(const double, double[], double[], int, const double[], const double[], const double[]);
typedef void (*external_func_keq_t)(int, const double[], double[], double[], const double);

} // namespace zerork

#endif

