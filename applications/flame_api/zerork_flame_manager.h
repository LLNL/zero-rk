#ifndef ZERORK_REACTOR_MANAGER_H
#define ZERORK_REACTOR_MANAGER_H


#include <string>
#include <vector>
#include <memory>
#include <fstream>

#include "zerork/mechanism.h"
#include <reactor/const_pressure_reactor.h>
#include <transport/mass_transport_factory.h>

#include "zerork_flame_api.h" //zerork_flame_status_t

#include "optionable.h"


class ZeroRKFlameManager : public Optionable
{
 public:
  ZeroRKFlameManager();
  virtual ~ZeroRKFlameManager() {};

  zerork_flame_status_t ReadOptionsFile(const std::string& options_filename);
  zerork_flame_status_t LoadMechanism();

//  zerork_flame_status_t SetCallbackFunction(zerork_callback_fn fn, void* cb_fn_data);

  zerork_flame_status_t FinishInit();
  zerork_flame_status_t Solve(int num_grid_points, const double* grid_points,
                              double P, double* flame_speed,
                              double* T, double* mass_fractions);

 private:
  std::shared_ptr<ConstPressureReactor> reactor_;
  std::shared_ptr<transport::MassTransportInterface> transport_;

  bool tried_init_;
  zerork_flame_status_t init_status_;
};

#endif
