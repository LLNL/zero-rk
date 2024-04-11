
//System libraries
#include <memory> //unique_ptr

//ZERORK
#include "zerork_cfd_plugin.h" //for zerork_handle, zerork_field_t
#include "zerork_reactor_manager.h"

#ifdef ZERORK_TRAP_FE
#include <fenv.h> //for fpe trapping
#endif

//#ifdef ZERORK_GPU
//#include "cuda_runtime_api.h"
//static void cuda_exit() {
//  cudaDeviceReset();
//}
//#endif

namespace { //anonymous
zerork_status_t OptionableStatusToZerorRKStatus(optionable_status_t status) {
  if(status == OPTIONABLE_STATUS_SUCCESS) return ZERORK_STATUS_SUCCESS;
  if(status == OPTIONABLE_STATUS_OPTION_NOT_DEFINED) return ZERORK_STATUS_OPTION_NOT_DEFINED;
  return ZERORK_STATUS_UNKNOWN_ERROR;
}
}

struct zerork_handle_impl {
  std::unique_ptr<ZeroRKReactorManagerBase> r;
};

extern "C"
zerork_handle zerork_reactor_init()
{
#ifdef ZERORK_TRAP_FE
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW );
#endif
//#ifdef ZERORK_GPU
//  atexit(cuda_exit);
//#endif

  zerork_handle handle = new zerork_handle_impl();
  handle->r = std::make_unique<ZeroRKReactorManager>();
  return handle;
}

extern "C"
zerork_status_t zerork_reactor_read_options_file(const char* options_filename,
                                     zerork_handle handle)
{
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  if(options_filename == nullptr) return ZERORK_STATUS_INVALID_FILENAME;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string options_filename_str(options_filename);
  return zrm->ReadOptionsFile(options_filename_str);
}

extern "C"
zerork_status_t zerork_reactor_set_mechanism_files(const char* mech_filename,
                                       const char* therm_filename,
                                       zerork_handle handle)
{
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  if(mech_filename == nullptr) return ZERORK_STATUS_INVALID_FILENAME;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = "mech_filename";
  std::string option_value_str = mech_filename;
  zerork_status_t flag = OptionableStatusToZerorRKStatus(zrm->SetStringOption(option_name, option_value_str));
  if(flag != ZERORK_STATUS_SUCCESS) return flag;
  if(therm_filename != nullptr) {
    option_name = "therm_filename";
    option_value_str = therm_filename;
    flag = OptionableStatusToZerorRKStatus(zrm->SetStringOption(option_name, option_value_str));
  }
  return flag;
}

zerork_status_t zerork_reactor_load_mechanism(zerork_handle handle)
{
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  return zrm->LoadMechanism();
}

extern "C"
zerork_status_t zerork_reactor_set_int_option(const char* option_name_chr,
                                  int option_value,
                                  zerork_handle handle)
{
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = option_name_chr;
  return OptionableStatusToZerorRKStatus(zrm->SetIntOption(option_name, option_value));
}

extern "C"
zerork_status_t zerork_reactor_set_double_option(const char* option_name_chr,
                                     double option_value,
                                     zerork_handle handle)
{
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = option_name_chr;
  return OptionableStatusToZerorRKStatus(zrm->SetDoubleOption(option_name, option_value));
}

extern "C"
zerork_status_t zerork_reactor_set_string_option(const char* option_name_chr,
                                     const char* option_value,
                                     zerork_handle handle)
{
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = option_name_chr;
  std::string option_value_str = option_value;
  return OptionableStatusToZerorRKStatus(zrm->SetStringOption(option_name, option_value_str));
}

extern "C"
zerork_status_t zerork_reactor_get_int_option(const char* option_name_chr,
                                  int* option_value,
                                  zerork_handle handle)
{
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = option_name_chr;
  return OptionableStatusToZerorRKStatus(zrm->GetIntOption(option_name, option_value));
}

extern "C"
zerork_status_t zerork_reactor_get_double_option(const char* option_name_chr,
                                     double* option_value,
                                     zerork_handle handle)
{
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = option_name_chr;
  return OptionableStatusToZerorRKStatus(zrm->GetDoubleOption(option_name, option_value));
}

extern "C"
zerork_status_t zerork_reactor_get_string_option(const char* option_name_chr,
                                     char** option_value,
                                     int string_length,
                                     zerork_handle handle)
{
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  std::string option_name = option_name_chr;
  std::string option_value_str;
  zerork_status_t flag = OptionableStatusToZerorRKStatus(zrm->GetStringOption(option_name, &option_value_str));
  if(flag == ZERORK_STATUS_SUCCESS) {
    strncpy(*option_value, option_value_str.c_str(), string_length);
  }
  return flag;
}

extern "C"
zerork_status_t zerork_reactor_solve(const int n_cycle, const double time,
                         const double dt, const int n_reactors,
                         double *T, double *P, double *mf,
                         zerork_handle handle)
{
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  //FinishInit may be first call to parse mechanism
  zerork_status_t flag = zrm->FinishInit();
  if(flag != ZERORK_STATUS_SUCCESS) return flag;
  zrm->SetInputVariables(n_cycle, time, dt, n_reactors, T, P, mf);
  zrm->LoadBalance();
  flag = zrm->SolveReactors();
  zrm->RedistributeResults();
  zrm->PostSolve();
  return flag;
}

extern "C"
zerork_status_t zerork_reactor_set_aux_field_pointer(zerork_field_t ft,
                                         double * field_pointer,
                                         zerork_handle handle) {
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  if(field_pointer == nullptr) return ZERORK_STATUS_INVALID_POINTER;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  return zrm->SetAuxFieldPointer(ft, field_pointer);
}

extern "C"
zerork_status_t zerork_reactor_set_callback_fn(zerork_callback_fn fn,
                                   void* fn_data,
                                   zerork_handle handle) {
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  if(fn == nullptr) return ZERORK_STATUS_INVALID_POINTER;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  return zrm->SetCallbackFunction(fn, fn_data);
}

extern "C"
zerork_status_t zerork_reactor_set_reactor_ids(int * reactor_ids,
                                   zerork_handle handle) {
  if(handle == nullptr) return ZERORK_STATUS_INVALID_HANDLE;
  if(reactor_ids == nullptr) return ZERORK_STATUS_INVALID_POINTER;
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  return zrm->SetReactorIDs(reactor_ids);
}

extern "C"
zerork_status_t zerork_reactor_free(zerork_handle handle) {
  if(handle == nullptr) {
      return ZERORK_STATUS_INVALID_HANDLE;
  } else {
      handle->r.reset(nullptr);
      delete handle;
      return ZERORK_STATUS_SUCCESS;
  }
}

