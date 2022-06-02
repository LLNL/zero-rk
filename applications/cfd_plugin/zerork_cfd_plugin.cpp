
//System libraries
#include <memory> //unique_ptr

//ZERORK
#include "zerork_cfd_plugin.h" //for zerork_handle, zerork_field_type
#include "zerork_reactor_manager.h"

#ifdef ZERORK_TRAP_FE
#include <fenv.h> //for fpe trapping
#endif

struct zerork_handle_impl {
  std::unique_ptr<ZeroRKReactorManagerBase> r;
};

extern "C"
zerork_handle zerork_reactor_init(const char* input_filename,
                        const char* mech_filename,
                        const char* therm_filename)
{
#ifdef ZERORK_TRAP_FE
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW );
#endif

  zerork_handle handle = new zerork_handle_impl();
  handle->r = std::make_unique<ZeroRKReactorManager>(input_filename, mech_filename, therm_filename);
  return handle;
}

extern "C"
int zerork_reactor_set_int_option(const char* option_name_chr,
                                  int option_value,
                                  zerork_handle handle)
{
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = option_name_chr;
  return zrm->SetIntOption(option_name, option_value);
}

extern "C"
int zerork_reactor_set_double_option(const char* option_name_chr,
                                     double option_value,
                                     zerork_handle handle)
{
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = option_name_chr;
  return zrm->SetDoubleOption(option_name, option_value);
}

extern "C"
int zerork_reactor_set_string_option(const char* option_name_chr,
                                     const char* option_value,
                                     zerork_handle handle)
{
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = option_name_chr;
  std::string option_value_str = option_value;
  return zrm->SetStringOption(option_name, option_value_str);
}

extern "C"
int zerork_reactor_get_int_option(const char* option_name_chr,
                                  int* option_value,
                                  zerork_handle handle)
{
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = option_name_chr;
  return zrm->GetIntOption(option_name, option_value);
}

extern "C"
int zerork_reactor_get_double_option(const char* option_name_chr,
                                     double* option_value,
                                     zerork_handle handle)
{
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = option_name_chr;
  return zrm->GetDoubleOption(option_name, option_value);
}

extern "C"
int zerork_reactor_get_string_option(const char* option_name_chr,
                                     char** option_value,
                                     int string_length,
                                     zerork_handle handle)
{
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = option_name_chr;
  std::string option_value_str;
  int flag = zrm->GetStringOption(option_name, &option_value_str);
  if(flag == 0) {
    strncpy(*option_value, option_value_str.c_str(), string_length);
  }
  return flag;
}

extern "C"
int zerork_reactor_solve(const int n_cycle, const double time,
                         const double dt, const int n_reactors,
                         double *T, double *P, double *mf,
                         zerork_handle handle)
{
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  zrm->FinishInit();
  zrm->SetInputVariables(n_cycle, time, dt, n_reactors, T, P, mf);
  zrm->LoadBalance();
  zrm->SolveReactors();
  zrm->RedistributeResults();
  zrm->PostSolve();
  return 0;
}

extern "C"
int zerork_reactor_set_aux_field_pointer(zerork_field_type ft,
                                         double * field_pointer,
                                         zerork_handle handle) {
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  zrm->SetAuxFieldPointer(ft, field_pointer);
  return 0;
}

extern "C"
int zerork_reactor_set_callback_fn(zerork_callback_fn fn,
                                   void* fn_data,
                                   zerork_handle handle) {
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  zrm->SetCallbackFunction(fn, fn_data);
  return 0;
}

extern "C"
int zerork_reactor_set_reactor_ids(int * reactor_ids,
                                   zerork_handle handle) {
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  zrm->SetReactorIDs(reactor_ids);
  return 0;
}

extern "C"
int zerork_reactor_free(zerork_handle handle) {
  if(handle == nullptr) {
      return 1;
  } else {
      handle->r.reset(nullptr);
      delete handle;
      return 0;
  }
}

