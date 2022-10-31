#ifndef FLAME_PARAMS_H_
#define FLAME_PARAMS_H_

#include <string>
#include <vector>

#include <file_utilities.h>

#ifdef ZERORK_MPI
#include <mpi.h>
#endif

#include <zerork/mechanism.h>
#include <reactor/const_pressure_reactor.h>
#include <transport/mass_transport_factory.h>

#include "sparse_matrix.h"

#include "UnsteadyFlameIFP.h"

class FlameParams
{
 public:
  explicit FlameParams(const std::string &input_name);
  ~FlameParams();

  // MPI
#ifdef ZERORK_MPI
  MPI_Comm comm_;
#endif
  int my_pe_,npes_;
  int num_procs_;
  int num_local_points_;
  int num_points_;

  zerork::mechanism *mechanism_; // TODO: avoid using a separate mechanism

  UnsteadyFlameIFP *parser_;
  ConstPressureReactor *reactor_;
  transport::MassTransportInterface *trans_;
  zerork::utilities::Logger *logger_;

  std::string input_name_;

  std::vector<int> fuel_species_id_;
  std::vector<int> oxidizer_species_id_;
  std::vector<int> full_species_id_;
  std::vector<int> egr_species_id_;

  std::vector<double> inlet_mass_fractions_;
  std::vector<double> initial_mass_fractions_;
  std::vector<double> z_;
  std::vector<double> zm_;
  std::vector<double> dz_;
  std::vector<double> dzm_;
  std::vector<double> dz_local_, dzm_local_, inv_dz_local_, inv_dzm_local_;
  std::vector<double> wall_temperature_;
  std::vector<double> mass_flux_;        // [kg/m^2/s] SetGrid()

  std::vector<double> y_ext_;
  std::vector<double> mass_flux_ext_;

  double mass_flux_inlet_;

  double diameter_;         // [m]
  double nusselt_;          // [-]        SetWallProperties()
  double inlet_temperature_;
  double inlet_molecular_mass_;
  double inlet_relative_volume_;
  double ref_temperature_;
  double flame_speed_;
  double flame_thickness_;
  double flame_thickness_alpha_;
  double max_temperature_;
  double mass_change_;
  double continuity_error_;

  double length_;

  std::vector<double> step_limiter_;

  // data to evaluate explicit time step limits, evaluated by the CVode
  // right hand side function
  double max_velocity_;          // [m/s]
  double max_thermal_diffusivity_;  // [m^2/s]

  double min_velocity_;          // [m/s]

  // data structures setup by SetMemory()
  transport::MassTransportInput transport_input_;
  std::vector<double> inv_molecular_mass_;     // size = num_species
  std::vector<double> species_specific_heats_; // size = num_species
  std::vector<double> species_mass_flux_;      // size = num_species
  std::vector<double> species_lewis_numbers_;      // size = num_species*num_points
  std::vector<double> thermal_conductivity_;   // size = num_points+1
  std::vector<double> mixture_specific_heat_;  // size = num_points
  std::vector<double> mixture_specific_heat_mid_;  // size = num_points+1

  int convective_scheme_type_; // 0 First order upwind
                               // 1 Second order upwind
                               // 2 Second order centered

  int simulation_type_; //0 for unsteady flame with extinction/ignition, 1 for steady flame

  // single reactor Jacobian structure (needed for sparse preconditioner
  // and banded solvers)

  bool store_jacobian_;
  bool valid_jacobian_structure_;
  std::vector<SparseMatrix *> sparse_matrix_;
  std::vector<int>     row_id_;
  std::vector<int>     column_sum_;
  std::vector<int>     diagonal_id_;
  std::vector<double>  reactor_jacobian_;
  std::vector<double>  saved_jacobian_;

 private:
  void SetInlet();
  void SetInitialComposition();
  void SetGrid();
  void SetWallProperties();
  void SetMemory();
};


#endif
