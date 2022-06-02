#ifndef ZERORK_CVREACTOR_LIB_H
#define ZERORK_CVREACTOR_LIB_H


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "zerork_cfd_plugin_exports.h"

typedef enum _zerork_field_type {
  ZERORK_FIELD_DPDT,
  ZERORK_FIELD_COST,
  ZERORK_FIELD_Y_SRC,
  ZERORK_FIELD_E_SRC,
  ZERORK_FIELD_IGNITION_TIME
} zerork_field_type;

typedef struct zerork_handle_impl* zerork_handle;
zerork_handle ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_init(const char* input_filename,
                        const char* mech_file,
                        const char* therm_file);

typedef int (*zerork_callback_fn)(int reactor_id, int steps, double time, double dt,
                                  const double* y, const double* dy,
                                  void* cb_fn_data);

int ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_solve(const int n_cycle, const double time,
                         const double dt, const int n_reactors,
                         double *T, double *P,
                         double *mf,
                         zerork_handle handle);

int ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_aux_field_pointer(zerork_field_type ft, double * field_pointer, zerork_handle handle);

int ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_int_option(const char* option_name_chr,
                                  int option_value,
                                  zerork_handle handle);

int ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_double_option(const char* option_name_chr,
                                     double option_value,
                                     zerork_handle handle);

int ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_string_option(const char* option_name_chr,
                                     const char* option_value,
                                     zerork_handle handle);

int ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_get_int_option(const char* option_name_chr,
                                  int* option_value,
                                  zerork_handle handle);

int ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_get_double_option(const char* option_name_chr,
                                     double* option_value,
                                     zerork_handle handle);

int ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_get_string_option(const char* option_name_chr,
                                     char** option_value,
                                     int string_length,
                                     zerork_handle handle);

int ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_callback_fn(zerork_callback_fn fn,
                                   void* cb_fn_data,
                                   zerork_handle handle);

int ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_reactor_ids(int* reactor_ids, zerork_handle handle);

int ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_free(zerork_handle handle); 


#ifdef __cplusplus
}
#endif


#endif
