#ifndef ZERORK_CVREACTOR_LIB_H
#define ZERORK_CVREACTOR_LIB_H


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "zerork_cfd_plugin_exports.h"

typedef enum _zerork_field_t {
  ZERORK_FIELD_DPDT,
  ZERORK_FIELD_COST,
  ZERORK_FIELD_GPU,
  ZERORK_FIELD_Y_SRC,
  ZERORK_FIELD_E_SRC,
  ZERORK_FIELD_IGNITION_TIME,
  ZERORK_FIELD_TEMPERATURE_DELTA
} zerork_field_t;

typedef enum _zerork_status_t {
  ZERORK_STATUS_SUCCESS = 0,
  ZERORK_STATUS_INVALID_HANDLE,
  ZERORK_STATUS_INVALID_FILENAME,
  ZERORK_STATUS_INVALID_FIELD_NAME,
  ZERORK_STATUS_INVALID_POINTER,
  ZERORK_STATUS_FAILED_MECHANISM_PARSE,
  ZERORK_STATUS_FAILED_OPTIONS_PARSE,
  ZERORK_STATUS_FAILED_SOLVE,
  ZERORK_STATUS_OPTION_NOT_DEFINED,
  ZERORK_STATUS_UNKNOWN_ERROR = 99,
} zerork_status_t;

typedef struct zerork_handle_impl* zerork_handle;
zerork_handle ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_init();

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_read_options_file(const char* options_filename, zerork_handle handle);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_mechanism_files(const char* mech_file, const char* therm_file,
                                                                             zerork_handle handle);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_load_mechanism(zerork_handle handle);

typedef int (*zerork_callback_fn)(int reactor_id, int steps, double time, double dt,
                                  const double* y, const double* dy,
                                  void* cb_fn_data);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_solve(const int n_cycle, const double time,
                         const double dt, const int n_reactors,
                         double *T, double *P,
                         double *mf,
                         zerork_handle handle);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_aux_field_pointer(zerork_field_t ft, double * field_pointer, zerork_handle handle);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_int_option(const char* option_name_chr,
                                  int option_value,
                                  zerork_handle handle);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_double_option(const char* option_name_chr,
                                     double option_value,
                                     zerork_handle handle);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_string_option(const char* option_name_chr,
                                     const char* option_value,
                                     zerork_handle handle);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_get_int_option(const char* option_name_chr,
                                  int* option_value,
                                  zerork_handle handle);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_get_double_option(const char* option_name_chr,
                                     double* option_value,
                                     zerork_handle handle);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_get_string_option(const char* option_name_chr,
                                     char** option_value,
                                     int string_length,
                                     zerork_handle handle);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_callback_fn(zerork_callback_fn fn,
                                   void* cb_fn_data,
                                   zerork_handle handle);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_set_reactor_ids(int* reactor_ids, zerork_handle handle);

zerork_status_t ZERORK_CFD_PLUGIN_EXPORTS zerork_reactor_free(zerork_handle handle);


#ifdef __cplusplus
}
#endif


#endif
