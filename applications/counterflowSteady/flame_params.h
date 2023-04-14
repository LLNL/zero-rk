#ifndef FLAME_PARAMS_H_
#define FLAME_PARAMS_H_

#include <string>
#include <vector>

#include <file_utilities.h>

#include <mpi.h>

#include <zerork/mechanism.h>
#include <reactor/counterflow_reactor.h>
#include <transport/mass_transport_factory.h>

#include "sparse_matrix.h"
#include "sparse_matrix_dist.h"

#include "SteadyFlameIFP.h"

#include <kinsol/kinsol.h>
//#include <nvector/nvector_parallel.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_math.h>


class FlameParams
{
 public:
  explicit FlameParams(const std::string &input_name);
  ~FlameParams();

  // MPI
  MPI_Comm comm_;
  MPI_Comm inter_comm_;
  int comm_rank_,num_comms_;
  int my_pe_,npes_;
  int num_procs_;
  int num_local_points_;
  int num_points_;
  int nover_;

  zerork::mechanism *mechanism_;

  SteadyFlameIFP *parser_;
  CounterflowReactor *reactor_;
  transport::MassTransportInterface *trans_;
  zerork::utilities::Logger *logger_;

  std::string input_name_;

  std::vector<int> fuel_species_id_;
  std::vector<int> oxidizer_species_id_;
  std::vector<int> full_species_id_;

  std::vector<double> fuel_mass_fractions_;
  std::vector<double> oxidizer_mass_fractions_;
  std::vector<double> stoichiometric_mass_fractions_;
  std::vector<double> inlet_mass_fractions_;
  std::vector<double> initial_mass_fractions_;
  std::vector<double> z_;
  std::vector<double> zm_;
  std::vector<double> dz_;
  std::vector<double> dzm_;
  std::vector<double> dz_local_, dzm_local_, inv_dz_local_, inv_dzm_local_;

  std::vector<double> y_ext_;
  std::vector<double> mass_flux_;
  std::vector<double> mass_flux_ext_;
  std::vector<double> rel_vol_;
  std::vector<double> rel_vol_ext_;
  std::vector<double> rhs_ext_;
  std::vector<double> rhsConv_;

  std::vector<double> y_old_;
  double dt_;

  double velocity_fuel_;
  double velocity_oxidizer_;

  double mass_flux_fuel_;
  double mass_flux_oxidizer_;

  double P_left_;
  double P_right_;

  double G_right_;

  double fuel_temperature_;
  double fuel_molecular_mass_;
  double fuel_relative_volume_;
  double oxidizer_temperature_;
  double oxidizer_molecular_mass_;
  double oxidizer_relative_volume_;
  double stoichiometric_molecular_mass_;
  double inlet_molecular_mass_;
  double inlet_relative_volume_;

  double ref_temperature_;
  double ref_momentum_;
  double flame_speed_;
  double flame_thickness_;
  double max_temperature_;

  double mass_change_;
  double continuity_error_;

  double stoichiometric_mixture_fraction_;

  double strain_rate_;
  double stagnation_plane_;
  int jContBC_;
  double length_;

  std::vector<double> step_limiter_;

  // data to evaluate explicit time step limits, evaluated by the CVode
  // right hand side function
  double max_velocity_;          // [m/s]
  double max_thermal_diffusivity_;  // [m^2/s]

  double min_velocity_;          // [m/s]

  double max_svf_;

  // data structures setup by SetMemory()
  transport::MassTransportInput transport_input_;
  std::vector<double> inv_molecular_mass_;     // size = num_species
  std::vector<double> species_specific_heats_; // size = num_species
  std::vector<double> species_mass_flux_;      // size = num_species
  std::vector<double> species_lewis_numbers_;      // size = num_species*num_points
  std::vector<double> thermal_conductivity_;   // size = num_points+1
  std::vector<double> mixture_viscosity_;   // size = num_points+1
  std::vector<double> mixture_specific_heat_;  // size = num_points
  std::vector<double> mixture_specific_heat_mid_;  // size = num_points+1
  std::vector<double> molecular_mass_mix_mid_;  // size = num_points+1

  //
  std::vector<double> fixed_temperature_; // size = num_points;

  int convective_scheme_type_; // 0 First order upwind
                               // 1 Second order upwind
                               // 2 Second order centered

  int simulation_type_; //0 for planar, 1 for axisymmetric
  int flame_type_; //0 for diffusion, 1 for premixed twin (R-to-R), 2 for premixed single (R-to-P)

  // single reactor Jacobian structure (needed for sparse preconditioner
  // and banded solvers)
  int integrator_type_;

  int num_off_diagonals_;

  bool store_jacobian_;
  bool valid_jacobian_structure_;

  bool pseudo_unsteady_;

  // For SuperLU serial
  bool superlu_serial_;
  SparseMatrix *sparse_matrix_;

  // For SuperLU_DIST
  SparseMatrix_dist *sparse_matrix_dist_;
  std::vector<double> reactor_jacobian_dist_;
  std::vector<double> saved_jacobian_dist_;
  std::vector<int>    col_id_;
  std::vector<int>    row_sum_;
  std::vector<int>    dense_to_sparse_;
  std::vector<int>    dense_to_sparse_offdiag_;
  int num_nonzeros_loc_;

  // For SuperLU serial for local chemistry
  std::vector<SparseMatrix *> sparse_matrix_chem_;
  std::vector<int>     row_id_chem_;
  std::vector<int>     column_id_chem_;
  std::vector<int>     column_sum_chem_;
  std::vector<int>     diagonal_id_chem_;
  std::vector<double>  reactor_jacobian_chem_;
  std::vector<double>  saved_jacobian_chem_;

  // For LAPACK serial banded solve
  int num_states_local_;
  int num_states_per_proc_;
  std::vector<double> banded_jacobian_;
  std::vector<double> banded_jacobian2_;
  std::vector<double> banded_jacobian_serial_;
  std::vector<int> pivots_serial_;

 private:
  void SetInlet();
  void SetGrid();
  void SetMemory();
};


#endif
