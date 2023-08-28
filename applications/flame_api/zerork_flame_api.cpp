
#include <memory> //unique_ptr


#ifndef _WIN32
#ifndef NDEBUG
#define ZERORK_TRAP_FE
#include <fenv.h> //for fpe trapping
#endif
#endif

#include "zerork_flame_api.h"
#include "zerork_flame_manager.h"

namespace { //anonymous
zerork_flame_status_t OptionableStatusToZerorRKStatus(optionable_status_t status) {
  if(status == OPTIONABLE_STATUS_SUCCESS) return ZERORK_FLAME_STATUS_SUCCESS;
  if(status == OPTIONABLE_STATUS_OPTION_NOT_DEFINED) return ZERORK_FLAME_STATUS_OPTION_NOT_DEFINED;
  return ZERORK_FLAME_STATUS_UNKNOWN_ERROR;
}
}

struct zerork_flame_handle_impl {
  std::unique_ptr<ZeroRKFlameManager> r;
};

extern "C"
zerork_flame_handle zerork_flame_init()
{
#ifdef ZERORK_TRAP_FE
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW );
#endif

  zerork_flame_handle handle = new zerork_flame_handle_impl();
  handle->r = std::make_unique<ZeroRKFlameManager>();
  return handle;
}

extern "C"
zerork_flame_status_t zerork_flame_read_options_file(const char* options_filename,
                                     zerork_flame_handle handle)
{
  if(handle == nullptr) return ZERORK_FLAME_STATUS_INVALID_HANDLE;
  if(options_filename == nullptr) return ZERORK_FLAME_STATUS_INVALID_FILENAME;
  ZeroRKFlameManager* zfm = handle->r.get();
  std::string options_filename_str(options_filename);
  return zfm->ReadOptionsFile(options_filename_str);
}


extern "C"
zerork_flame_status_t zerork_flame_set_input_files(const char* mech_filename,
                                                   const char* therm_filename,
                                                   const char* trans_filename,
                                                   zerork_flame_handle handle)
{
  if(handle == nullptr) return ZERORK_FLAME_STATUS_INVALID_HANDLE;
  if(mech_filename == nullptr) return ZERORK_FLAME_STATUS_INVALID_FILENAME;
  if(trans_filename == nullptr) return ZERORK_FLAME_STATUS_INVALID_FILENAME;
  ZeroRKFlameManager* zfm = handle->r.get();
  std::string option_name = "mechanism_filename";
  std::string option_value_str = mech_filename;
  zerork_flame_status_t flag = OptionableStatusToZerorRKStatus(zfm->SetStringOption(option_name, option_value_str));
  if(flag != ZERORK_FLAME_STATUS_SUCCESS) return flag;
  if(therm_filename != nullptr) {
    option_name = "therm_filename";
    option_value_str = therm_filename;
    flag = OptionableStatusToZerorRKStatus(zfm->SetStringOption(option_name, option_value_str));
    if(flag != ZERORK_FLAME_STATUS_SUCCESS) return flag;
  }
  option_name = "transport_filename";
  option_value_str = trans_filename;
  flag = OptionableStatusToZerorRKStatus(zfm->SetStringOption(option_name, option_value_str));
  if(flag != ZERORK_FLAME_STATUS_SUCCESS) return flag;
  return flag;
}

extern "C"
zerork_flame_status_t zerork_flame_load_mechanism(zerork_flame_handle handle)
{
  if(handle == nullptr) return ZERORK_FLAME_STATUS_INVALID_HANDLE;
  ZeroRKFlameManager* zfm = handle->r.get();
  return zfm->LoadMechanism();
}

extern "C"
zerork_flame_status_t zerork_flame_solve(int num_grid_points, const double* grid_points,
                         double P, double* flame_speed, double* T, double* mass_fractions,
                         zerork_flame_handle handle)
{
  if(handle == nullptr) return ZERORK_FLAME_STATUS_INVALID_HANDLE;
  ZeroRKFlameManager* zfm = handle->r.get();
  //FinishInit may be first call to parse mechanism
  zerork_flame_status_t flag = zfm->FinishInit();
  if(flag != ZERORK_FLAME_STATUS_SUCCESS) return flag;
  flag = zfm->Solve(num_grid_points, grid_points, P, flame_speed, T, mass_fractions);
  return flag;
}

extern "C"
zerork_flame_status_t zerork_flame_set_int_option(const char* option_name_chr,
                                  int option_value,
                                  zerork_flame_handle handle)
{
  if(handle == nullptr) return ZERORK_FLAME_STATUS_INVALID_HANDLE;
  ZeroRKFlameManager* zfm = handle->r.get();
  std::string option_name = option_name_chr;
  return OptionableStatusToZerorRKStatus(zfm->SetIntOption(option_name, option_value));
}

extern "C"
zerork_flame_status_t zerork_flame_set_double_option(const char* option_name_chr,
                                     double option_value,
                                     zerork_flame_handle handle)
{
  if(handle == nullptr) return ZERORK_FLAME_STATUS_INVALID_HANDLE;
  ZeroRKFlameManager* zfm = handle->r.get();
  std::string option_name = option_name_chr;
  return OptionableStatusToZerorRKStatus(zfm->SetDoubleOption(option_name, option_value));
}

extern "C"
zerork_flame_status_t zerork_flame_set_string_option(const char* option_name_chr,
                                     const char* option_value,
                                     zerork_flame_handle handle)
{
  if(handle == nullptr) return ZERORK_FLAME_STATUS_INVALID_HANDLE;
  ZeroRKFlameManager* zfm = handle->r.get();
  std::string option_name = option_name_chr;
  std::string option_value_str = option_value;
  return OptionableStatusToZerorRKStatus(zfm->SetStringOption(option_name, option_value_str));
}

extern "C"
zerork_flame_status_t zerork_flame_get_int_option(const char* option_name_chr,
                                  int* option_value,
                                  zerork_flame_handle handle)
{
  if(handle == nullptr) return ZERORK_FLAME_STATUS_INVALID_HANDLE;
  ZeroRKFlameManager* zfm = handle->r.get();
  std::string option_name = option_name_chr;
  return OptionableStatusToZerorRKStatus(zfm->GetIntOption(option_name, option_value));
}

extern "C"
zerork_flame_status_t zerork_flame_get_double_option(const char* option_name_chr,
                                     double* option_value,
                                     zerork_flame_handle handle)
{
  if(handle == nullptr) return ZERORK_FLAME_STATUS_INVALID_HANDLE;
  ZeroRKFlameManager* zfm = handle->r.get();
  std::string option_name = option_name_chr;
  return OptionableStatusToZerorRKStatus(zfm->GetDoubleOption(option_name, option_value));
}

extern "C"
zerork_flame_status_t zerork_flame_get_string_option(const char* option_name_chr,
                                     char** option_value,
                                     int string_length,
                                     zerork_flame_handle handle)
{
  if(handle == nullptr) return ZERORK_FLAME_STATUS_INVALID_HANDLE;
  ZeroRKFlameManager* zfm = handle->r.get();
  std::string option_name = option_name_chr;
  std::string option_value_str;
  zerork_flame_status_t flag = OptionableStatusToZerorRKStatus(zfm->GetStringOption(option_name, &option_value_str));
  if(flag == ZERORK_FLAME_STATUS_SUCCESS) {
    strncpy(*option_value, option_value_str.c_str(), string_length);
  }
  return flag;
}

extern "C"
zerork_flame_status_t zerork_flame_free(zerork_flame_handle handle) {
  if(handle == nullptr) {
      return ZERORK_FLAME_STATUS_INVALID_HANDLE;
  } else {
      handle->r.reset(nullptr);
      delete handle;
      return ZERORK_FLAME_STATUS_SUCCESS;
  }
}




