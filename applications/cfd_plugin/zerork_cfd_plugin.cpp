
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

namespace {
  std::vector<std::unique_ptr<zerork_handle_impl> > reactor_manager_handles;
}

extern "C"
zerork_handle zerork_reactor_init(const char* input_filename,
                        const char* mech_filename,
                        const char* therm_filename)
{
#ifdef ZERORK_TRAP_FE
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW );
#endif

  //TODO: Thread safety
  reactor_manager_handles.emplace_back(std::make_unique<zerork_handle_impl>());
  reactor_manager_handles.back()->r = std::make_unique<ZeroRKReactorManager>(input_filename, mech_filename, therm_filename);
  zerork_handle handle = reactor_manager_handles.back().get();
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
int zerork_reactor_solve(const int n_cycle, const double time,
                         const double dt, const int n_reactors,
                         double *T, double *P, double *mf,
                         zerork_handle handle)
{
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  zrm->SetInputVariables(n_cycle, time, dt, n_reactors, T, P, mf);
  zrm->LoadBalance();
  zrm->SolveReactors();
  zrm->RedistributeResults();
  zrm->PostSolve();
  return 0;
}

int zerork_reactor_set_aux_field_pointer(zerork_field_type ft,
                                         double * field_pointer,
                                         zerork_handle handle) {
  ZeroRKReactorManagerBase* zrm = handle->r.get();
  zrm->SetAuxFieldPointer(ft, field_pointer);
  return 0;
}

extern "C"
int zerork_reactor_free(zerork_handle handle) {
  if(handle == nullptr) {
      return 1;
  } else {
      handle->r.reset(nullptr);
      return 0;
  }
}

