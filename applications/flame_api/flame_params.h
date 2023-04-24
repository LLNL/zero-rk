#ifndef FLAME_PARAMS_H_
#define FLAME_PARAMS_H_

#ifdef ZERORK_MPI
#include <mpi.h>
#endif

#include <vector>

#include "optionable.h"

#include <zerork/mechanism.h>
#include <reactor/const_pressure_reactor.h>
#include <transport/mass_transport_factory.h>

#include "sparse_matrix.h"
#ifdef ZERORK_MPI
#include "sparse_matrix_dist.h"
#endif

class FlameParams : public Optionable
{
 public:
  explicit FlameParams(ConstPressureReactor* reactor, 
                       transport::MassTransportInterface* trans,
                       const std::vector<double>& grid,
                       double flame_speed,
                       const double* T,
                       const double* mass_fractions,
                       double pressure,
#ifdef ZERORK_MPI
                       MPI_Comm comm,
#endif
                       const Optionable options);
  ~FlameParams();

#ifdef ZERORK_MPI
  MPI_Comm comm_;
#endif
  int error_status_;
  int my_pe_,npes_;
  int num_procs_;
  int num_local_points_;
  int num_points_;
  int num_species_;
  int num_states_;
  int nover_;
  int num_kinsol_errors_;

  ConstPressureReactor *reactor_;
  transport::MassTransportInterface *transport_;
  zerork::mechanism *mechanism_;

  std::vector<double> inlet_mass_fractions_;
  std::vector<double> z_;
  std::vector<double> zm_;
  std::vector<double> dz_;
  std::vector<double> dzm_;
  std::vector<double> dz_local_, dzm_local_, inv_dz_local_, inv_dzm_local_;
  std::vector<double> rel_vol_;        // [kg/m^2/s] SetGrid()

  std::vector<double> y_ext_;
  std::vector<double> rhs_ext_;
  std::vector<double> rhsConv_;
  std::vector<double> rel_vol_ext_;

  std::vector<double> y_;
  std::vector<double> y_old_;
  double dt_;
  double pressure_;

  double mass_flux_;
  double mass_flux_inlet_;

  double inlet_temperature_;
  double outlet_temperature_;
  double inlet_relative_volume_;
  double reference_temperature_;
  double flame_speed_;
  double flame_thickness_;
  double flame_thickness_alpha_;
  double max_temperature_;
  double mass_change_;
  double continuity_error_;

  double temperature_fix_;
  int j_fix_;

  // data structures setup by SetMemory()
  transport::MassTransportInput transport_input_;
  //Containers for MassTransportInput 
  std::vector<double> mti_mf_;        // size = num_species_
  std::vector<double> mti_mf_grad_;   // size = num_species_
  std::vector<double> mti_temp_grad_; // size = num_dim
  std::vector<double> mti_pres_grad_; // size = num_dim

  std::vector<double> inv_molecular_mass_;     // size = num_species
  std::vector<double> species_specific_heats_; // size = num_species
  std::vector<double> species_mass_flux_;      // size = num_species
  std::vector<double> species_lewis_numbers_;      // size = num_species*num_points
  std::vector<double> thermal_conductivity_;   // size = num_points+1
  std::vector<double> mixture_specific_heat_;  // size = num_points
  std::vector<double> mixture_specific_heat_mid_;  // size = num_points+1
  std::vector<double> molecular_mass_mix_mid_;  // size = num_points+1

  int convective_scheme_type_; // 0 First order upwind
                               // 1 Second order upwind
                               // 2 Second order centered

  int num_off_diagonals_;
  int storage_mu_;

  // single reactor Jacobian structure (needed for sparse preconditioner
  // and banded solvers
  int integrator_type_; // 0 KINSOL BBD
                        // 2 ScaLAPACK banded matrix
                        // 3 SuperLU sparse matrix (global)
                        // 4 SuperLU sparse (local) + LAPACK tridiagonal
  bool store_jacobian_;
  bool valid_jacobian_structure_;

  bool pseudo_unsteady_;

  std::vector<double> step_limiter_;

  bool superlu_serial_;

  // For SuperLU serial
  SparseMatrix *sparse_matrix_;

  // For SuperLU_DIST
#ifdef ZERORK_MPI
  SparseMatrix_dist *sparse_matrix_dist_;
#endif
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

  void SetTfix();
  void GetTemperatureAndMassFractions(double* T, double* mass_fractions);
 private:
  void SetInitialCondition(double flame_speed, const double* T, const double* mass_fractions);
  void SetGrid();
  void SetMemory();
};


#endif
