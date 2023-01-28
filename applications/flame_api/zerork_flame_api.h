#ifndef ZERORK_FLAME_H
#define ZERORK_FLAME_H


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "zerork_flame_exports.h"

typedef enum _zerork_flame_status_t {
  ZERORK_FLAME_STATUS_SUCCESS = 0,
  ZERORK_FLAME_STATUS_INVALID_HANDLE,
  ZERORK_FLAME_STATUS_INVALID_FILENAME,
  ZERORK_FLAME_STATUS_INVALID_POINTER,
  ZERORK_FLAME_STATUS_FAILED_MECHANISM_PARSE,
  ZERORK_FLAME_STATUS_FAILED_TRANSPORT_PARSE,
  ZERORK_FLAME_STATUS_FAILED_OPTIONS_PARSE,
  ZERORK_FLAME_STATUS_FAILED_SOLVE,
  ZERORK_FLAME_STATUS_OPTION_NOT_DEFINED,
  ZERORK_FLAME_STATUS_UNKNOWN_ERROR,
} zerork_flame_status_t;

typedef struct zerork_flame_handle_impl* zerork_flame_handle;
zerork_flame_handle ZERORK_FLAME_EXPORTS zerork_flame_init();

zerork_flame_status_t ZERORK_FLAME_EXPORTS zerork_flame_read_options_file(const char* options_filename, zerork_flame_handle handle);

zerork_flame_status_t ZERORK_FLAME_EXPORTS zerork_flame_set_input_files(const char* mech_file, const char* therm_file, const char* trans_file,
                                                                        zerork_flame_handle handle);

zerork_flame_status_t ZERORK_FLAME_EXPORTS zerork_flame_load_mechanism(zerork_flame_handle handle);

//GT:  What information/variables do you want passed into the callback function?
//typedef int (*zerork_flame_callback_fn)(int steps, double time, double dt, const double flame_speed,
//                                        const double* T, const double* mf, void* cb_fn_data);


//GT: We can separate out the grid-info and pressure from this call, but decided this was simpler
zerork_flame_status_t ZERORK_FLAME_EXPORTS zerork_flame_solve(int num_grid_points,       //scalar (in)
                                                              const double *grid_points, //array (in)
                                                              double P,                  //scalar (in)
                                                              double *flame_speed,       //scalar (in/out)
                                                              double *T,                 //array (in/out) of temperatures at each grid point
                                                              double *mass_fractions,    //array (in/out) of mass fractions ([grid_idx*num_species + species_idx])
                                                              zerork_flame_handle handle);

zerork_flame_status_t ZERORK_FLAME_EXPORTS zerork_flame_set_int_option(const char* option_name_chr,
                                  int option_value,
                                  zerork_flame_handle handle);

zerork_flame_status_t ZERORK_FLAME_EXPORTS zerork_flame_set_double_option(const char* option_name_chr,
                                     double option_value,
                                     zerork_flame_handle handle);

//Transport model (const-Le vs. mixture-average) should be set here
zerork_flame_status_t ZERORK_FLAME_EXPORTS zerork_flame_set_string_option(const char* option_name_chr,
                                     const char* option_value,
                                     zerork_flame_handle handle);

zerork_flame_status_t ZERORK_FLAME_EXPORTS zerork_flame_get_int_option(const char* option_name_chr,
                                  int* option_value,
                                  zerork_flame_handle handle);

zerork_flame_status_t ZERORK_FLAME_EXPORTS zerork_flame_get_double_option(const char* option_name_chr,
                                     double* option_value,
                                     zerork_flame_handle handle);

zerork_flame_status_t ZERORK_FLAME_EXPORTS zerork_flame_get_string_option(const char* option_name_chr,
                                     char** option_value,
                                     int string_length,
                                     zerork_flame_handle handle);

//zerork_flame_status_t ZERORK_FLAME_EXPORTS zerork_flame_set_callback_fn(zerork_flame_callback_fn fn,
//                                                                        void* cb_fn_data,
//                                                                        zerork_flame_handle handle);

zerork_flame_status_t ZERORK_FLAME_EXPORTS zerork_flame_free(zerork_flame_handle handle);


#ifdef __cplusplus
}
#endif


#endif
