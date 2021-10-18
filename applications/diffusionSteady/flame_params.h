#ifndef FLAME_PARAMS_H_
#define FLAME_PARAMS_H_

#include <string>
#include <vector>

#include <file_utilities.h>

#include <mpi.h>

#include <zerork/mechanism.h>
#include <reactor/const_pressure_reactor.h>
#include <transport/mass_transport_factory.h>

#include "sparse_matrix.h"
#include "sparse_matrix_dist.h"

#include "UnsteadyFlameIFP.h"

#include <kinsol/kinsol.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_math.h>

#include "soot.h"

class FlameParams
{
 public:
  explicit FlameParams(const std::string &input_name, MPI_Comm &comm);
  ~FlameParams();

  //mpi
  MPI_Comm comm_;
  int my_pe_,npes_;
  int num_procs_;
  int num_local_points_;
  int num_points_;
  int nover_;

  UnsteadyFlameIFP *parser_;
  ConstPressureReactor *reactor_;
  transport::MassTransportInterface *trans_;
  zerork::utilities::Logger *logger_;

  zerork::mechanism *mechanism_; // TODO: avoid using a separate mechanism

  std::string input_name_;

  std::vector<int> fuel_species_id_;
  std::vector<int> oxidizer_species_id_;

  std::vector<double> fuel_mass_fractions_;
  std::vector<double> oxidizer_mass_fractions_;
  std::vector<double> stoichiometric_mass_fractions_;
  std::vector<double> initial_mass_fractions_;
  std::vector<double> z_;
  std::vector<double> zm_;
  std::vector<double> dz_;
  std::vector<double> dzm_;
  std::vector<double> dz_local_, dzm_local_, inv_dz_local_, inv_dzm_local_;
  std::vector<double> fixed_temperature_;
  bool fix_temperature_;
  std::vector<double> y_ext_;
  std::vector<double> rhs_ext_;

  std::vector<double> y_old_;
  double dt_;

  double fuel_temperature_;
  double fuel_molecular_mass_;
  double fuel_relative_volume_;
  double oxidizer_temperature_;
  double oxidizer_molecular_mass_;
  double oxidizer_relative_volume_;
  double stoichiometric_molecular_mass_;
  double ref_temperature_;
  double flame_thickness_;
  double max_temperature_;

  double stoichiometric_mixture_fraction_;
  double scalar_dissipation_rate_;

  std::vector<double> step_limiter_;

  double jacobian_constant_;
  int krylov_subspace_;

  // data to evaluate explicit time step limits, evaluated by the CVode
  // right hand side function
  double max_thermal_diffusivity_;  // [m^2/s]

  // data structures setup by SetMemory()
  transport::MassTransportInput transport_input_;
  std::vector<double> inv_molecular_mass_;     // size = num_species
  std::vector<double> species_specific_heats_; // size = num_species
  std::vector<double> species_mass_flux_;      // size = num_species
  std::vector<double> species_lewis_numbers_;      // size = num_species*(num_points+2)
  std::vector<double> thermal_conductivity_;   // size = num_points+2
  std::vector<double> mixture_specific_heat_;  // size = num_points+2

  std::vector<double> dissipation_rate_; //size = num_points+2
  std::vector<double> enthalpy_flux_sum_;
  std::vector<double> molecular_mass_;
  std::vector<double> mixture_molecular_mass_;
  std::vector<double> mixture_molecular_mass_grad_;
  std::vector<double> mixture_molecular_mass_laplacian_;
  std::vector<double> sum_mass_fraction_over_Lewis_;
  std::vector<double> sum_mass_fraction_over_Lewis_grad_;
  std::vector<double> sum_mass_fraction_over_Lewis_laplacian_;
  std::vector<double> convection_velocity_;
  std::vector<double> rho_dot_;

  // Smooth values for stability?
  std::vector<double> mixture_molecular_mass_smooth_;
  std::vector<double> sum_mass_fraction_over_Lewis_smooth_;

  int convective_scheme_type_; // 0 First order upwind
                               // 1 Second order upwind
                               // 2 Second order centered

  int simulation_type_; //0 for FREI, 1 for steady flame

  // single reactor Jacobian structure (needed for sparse preconditioner
  // and banded solvers
  int integrator_type_; // 0 KINSOL BBD
                        // 2 SuperLU sparse matrix
                        // 3 SuperLU sparse + LAPACK tridiagonal
  int num_off_diagonals_;

  bool store_jacobian_;
  bool valid_jacobian_structure_;

  bool pseudo_unsteady_;

  bool unity_Lewis_;
  bool full_equations_;
  bool soot_;
  bool sensitivity_analysis_;
  bool uncertainty_quantification_;

  // For SuperLU serial for local chemistry
  std::vector<SparseMatrix *> sparse_matrix_;
  std::vector<int>     row_id_;
  std::vector<int>     column_id_;
  std::vector<int>     column_sum_;
  std::vector<int>     diagonal_id_;
  std::vector<double>  reactor_jacobian_;
  std::vector<double>  saved_jacobian_;

  std::vector<double>  saved_jacobian_product_;

  // For SuperLU_DIST
  SparseMatrix_dist *sparse_matrix_dist_;
  std::vector<double> reactor_jacobian_dist_;
  std::vector<double> saved_jacobian_dist_;
  std::vector<int>    col_id_;
  std::vector<int>    row_sum_;
  std::vector<int>    dense_to_sparse_;
  std::vector<int>    dense_to_sparse_offdiag_;
  int num_nonzeros_loc_;

  // For LAPACK serial banded solve
  int num_states_local_;
  int num_states_per_proc_;
  std::vector<double> banded_jacobian_;
  std::vector<double> banded_jacobian2_;
  std::vector<double> banded_jacobian_serial_;
  std::vector<int> pivots_serial_;

  void SetInletBL();
  void SetInlet();

 private:
  void SetGrid();
  void SetFixedTProperties();
  void SetMemory();
};


#endif
