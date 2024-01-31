#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
#include <algorithm>
#include <utilities/string_utilities.h>
#include <utilities/math_utilities.h>
#include <utilities/file_utilities.h>

#include "flame_params.h"
#include "transport/flexible_transport.h"

FlameParams::FlameParams(ConstPressureReactor* reactor, 
                         transport::MassTransportInterface* trans,
                         const std::vector<double>& grid,
                         double flame_speed,
                         const double* T,
                         const double* mass_fractions,
                         const double pressure,
#ifdef ZERORK_MPI
                         MPI_Comm comm,
#endif
                         const Optionable options) :
  reactor_(reactor),
  transport_(trans),
  num_points_(grid.size()-1),
  z_(grid),
  num_states_(reactor->GetNumStates()),
  num_species_(reactor->GetNumSpecies()),
  num_kinsol_errors_(0),
  pressure_(pressure),
#ifdef ZERORK_MPI
  comm_(comm),
#endif
  Optionable(options),
  mechanism_(reactor->GetMechanism())
{
  error_status_ = 0;
  npes_ = 1;
  my_pe_ = 0;
#ifdef ZERORK_MPI
  MPI_Comm_size(comm_, &npes_);
  MPI_Comm_rank(comm_, &my_pe_);
#endif
  nover_ = 2;
  sparse_matrix_ = NULL;
#ifdef ZERORK_MPI
  sparse_matrix_dist_ = NULL;
#endif
  valid_jacobian_structure_ = true;

  integrator_type_ = int_options_["integrator_type"];
  if(integrator_type_ < 2 || integrator_type_ > 3) {
    error_status_ = 1;
  }
  convective_scheme_type_ = int_options_["convective_scheme_type"];
  if(convective_scheme_type_ < 0  || convective_scheme_type_ > 2) {
    error_status_ = 1;
  }
  store_jacobian_  = int_options_["store_jacobian"] != 0;
  reference_temperature_ =  double_options_["reference_temperature"];
  step_limiter_.assign(reactor_->GetNumSteps(), double_options_["step_limiter"]);
  pseudo_unsteady_ = int_options_["pseudo_unsteady"] != 0;
  dt_ = double_options_["pseudo_unsteady_dt"];
  if(dt_ <= 0.0) {
    error_status_ = 1;
  }

  if(error_status_ == 0) {
    SetInitialCondition(flame_speed,T,mass_fractions);
  }
  if(error_status_ == 0) {
    SetGrid();
  }
  if(error_status_ == 0) {
    SetTfix();
  }
  if(error_status_ == 0) {
    SetMemory();
  }
  if(error_status_ == 0) {
    if(string_options_["transport_model"] == "ConstantLewis") {
      bool need_lewis_update = false;
      int lewis_grid_point = -1;
      if(string_options_["constant_lewis_setting"] == "GridPoint") {
        need_lewis_update = true;
        lewis_grid_point = int_options_["constant_lewis_grid_point"];
        if(lewis_grid_point == -1) {
          lewis_grid_point = num_points_ - 1;
        }
      } else if(string_options_["constant_lewis_setting"] == "Tfix") {
        need_lewis_update = true;
        lewis_grid_point = j_fix_;
      }
      if(need_lewis_update) {
        lewis_grid_point = std::min(std::max(0,lewis_grid_point),num_points_-1);

#ifdef ZERORK_MPI
        int lewis_grid_point_rank = lewis_grid_point / num_local_points_;
        int lewis_grid_point_local = lewis_grid_point - (lewis_grid_point_rank*num_local_points_);
        std::vector<double> lewis_state(num_states_);
        if(my_pe_ == lewis_grid_point_rank) {
          for(int k=0; k<num_states_; ++k) {
            lewis_state[k] = y_[lewis_grid_point_local*num_states_+k];
          }
        }
        MPI_Bcast(&lewis_state[0], num_states_, MPI_DOUBLE, lewis_grid_point_rank, comm_);

        for(int k=0; k<num_species_; ++k) {
          transport_input_.mass_fraction_[k] = lewis_state[k];
          transport_input_.grad_mass_fraction_[k] = 0.0;
        }
        transport_input_.temperature_ = reference_temperature_*lewis_state[num_states_-1];
        transport_input_.grad_temperature_[0] = 0.0;
#else
        for(int k=0; k<num_species_; ++k) {
          transport_input_.mass_fraction_[k] = y_[lewis_grid_point*num_states_+k];
          transport_input_.grad_mass_fraction_[k] = 0.0;
        }
        transport_input_.temperature_ = reference_temperature_*y_[(lewis_grid_point+1)*num_states_-1];
        transport_input_.grad_temperature_[0] = 0.0;
#endif

        transport::FlexibleTransport* flexible_transport = dynamic_cast<transport::FlexibleTransport*>(transport_);
        flexible_transport->SetMixAvg(true);
        int transport_error = flexible_transport->GetSpeciesMassFlux(
				   transport_input_,
				   num_species_,
				   nullptr, //conductivity (computed internally)
				   nullptr, //specific heat (computed(internally)
				   &species_mass_flux_[0], //unused
				   &species_lewis_numbers_[0]); //unused
        flexible_transport->SetMixAvg(false);
        if(transport_error != transport::NO_ERROR) {
          error_status_ = 1;
        }
      }
    }
  }
}

FlameParams::~FlameParams()
{
  if(integrator_type_ == 2) {
    if (superlu_serial_) {
      if (sparse_matrix_ != NULL) {
        sparse_matrix_->SparseMatrixClean();
	delete sparse_matrix_;
      }
#ifdef ZERORK_MPI
    } else {
      if (sparse_matrix_dist_ != NULL) {
        sparse_matrix_dist_->SparseMatrixClean_dist();
	delete sparse_matrix_dist_;
      }
#endif
    }
  }
  if(integrator_type_ == 3) {
    for(size_t j=0; j<sparse_matrix_chem_.size(); ++j) {
      if(sparse_matrix_chem_[j] != NULL) {
        sparse_matrix_chem_[j]->SparseMatrixClean();
	delete sparse_matrix_chem_[j];
      }
    }
  }
}

void FlameParams::SetInitialCondition(double flame_speed, const double* T, const double* mass_fractions)
{
#ifdef ZERORK_MPI
  //Distribute inputs over ranks
  MPI_Bcast(&num_points_, 1, MPI_INT, 0, comm_);
  num_local_points_ = num_points_/npes_;
  int remainder = num_points_ % npes_;
  if(remainder != 0) {
    error_status_ = 1;
    return;
  }
  if(my_pe_ < remainder) {
    num_local_points_ += 1;
  }
  std::vector<int> scatter_counts(npes_,0);
  std::vector<int> scatter_counts_mf(npes_,0);
  std::vector<int> scatter_displs(npes_,0);
  std::vector<int> scatter_displs_mf(npes_,0);
  for(int j = 0; j < npes_; ++j) {
    scatter_counts[j] = num_points_/npes_;
    if(j<remainder) scatter_counts[j] += 1;
    scatter_counts_mf[j] = scatter_counts[j]*num_species_;
    if(j>0) {
      scatter_displs[j] = scatter_displs[j-1]+scatter_counts[j-1];
      scatter_displs_mf[j] = scatter_displs_mf[j-1]+scatter_counts_mf[j-1];
    }
  }

  std::vector<double> mass_fractions_local(num_local_points_*num_species_,0.0);
  std::vector<double> T_local(num_local_points_,0.0);
  MPI_Scatterv(&mass_fractions[num_species_], scatter_counts_mf.data(), scatter_displs_mf.data(), MPI_DOUBLE,
               mass_fractions_local.data(), scatter_counts_mf[my_pe_], MPI_DOUBLE, 0, comm_);
  MPI_Scatterv(&T[1], scatter_counts.data(), scatter_displs.data(), MPI_DOUBLE,
               T_local.data(), scatter_counts[my_pe_], MPI_DOUBLE, 0, comm_);
#else
  num_local_points_ = num_points_;
  const double* mass_fractions_local = &mass_fractions[num_species_];
  const double* T_local = &T[1];
#endif
  y_.resize(num_local_points_*num_states_);

  for(int j = 0; j < num_local_points_; ++j) {
    for(int k = 0; k < num_species_; ++k) {
      y_[j*num_states_ + k] = mass_fractions_local[j*num_species_ + k];
    }
    y_[j*num_states_ + num_species_ + 1] = T_local[j]/reference_temperature_;
  }

  if(my_pe_ == 0) {
    std::vector<double> inlet_mole_fractions(num_species_, 0.0);
    inlet_mass_fractions_.assign(num_species_,0.0);
    for(int j=0; j<num_species_; ++j) {
      inlet_mass_fractions_[j] = mass_fractions[j];
    }

    mechanism_->getXfromY(inlet_mass_fractions_.data(), inlet_mole_fractions.data());
    double inlet_molecular_mass = mechanism_->getMolWtMixFromX(inlet_mole_fractions.data());

    inlet_temperature_ = T[0]/reference_temperature_;
    inlet_relative_volume_ = reactor_->GetGasConstant()*inlet_temperature_*reference_temperature_/
                               (pressure_*inlet_molecular_mass);
    //S_L = y[num_species_]*params->inlet_relative_volume_
    mass_flux_ = flame_speed/inlet_relative_volume_;
  }
  if(my_pe_ == npes_-1) {
    outlet_temperature_ = y_[(num_local_points_-1)*num_states_ + num_species_+1];
  }
#ifdef ZERORK_MPI
  //N.B. inlet_temperature_ used in jfix Tfix, inlet_relative_volume_ currently only used on root, 
  //     but Bcast to avoid issues if it's needed elsewhere in the future.
  MPI_Bcast(&inlet_temperature_, 1, MPI_DOUBLE, 0, comm_);
  MPI_Bcast(&inlet_relative_volume_, 1, MPI_DOUBLE, 0, comm_);
  MPI_Bcast(&mass_flux_, 1, MPI_DOUBLE, 0, comm_);
  MPI_Bcast(&outlet_temperature_, 1, MPI_DOUBLE, npes_-1, comm_);
#endif

  for(int j = 0; j < num_local_points_; ++j) {
    y_[j*num_states_ + num_species_] = mass_flux_;
  }
} // void FlameParams::SetInitialCondition()

// Set the grid
void FlameParams::SetGrid()
{
  if(my_pe_ == 0) {
    z_.erase(z_.begin());
  }
#ifdef ZERORK_MPI
  z_.resize(num_points_);
  MPI_Bcast(z_.data(), num_points_, MPI_DOUBLE, 0, comm_);
#endif

  // Compute midpoints and spacings
  zm_.assign(num_points_, 0.0);
  dz_.assign(num_points_, 1.0);
  dzm_.assign(num_points_, 1.0);

  if(num_local_points_ < nover_) {
    error_status_ = 2;
    return;
  }

  dz_local_.assign( num_local_points_+(2*nover_), 0.0);
  dzm_local_.assign( num_local_points_+(2*nover_), 0.0);
  inv_dz_local_.assign( num_local_points_+(2*nover_), 0.0);
  inv_dzm_local_.assign( num_local_points_+(2*nover_), 0.0);

  // Compute midpoints and grid spacings
  for(int j=1; j<num_points_; ++j) {
    dz_[j] = z_[j]-z_[j-1];
    zm_[j] = 0.5*(z_[j]+z_[j-1]);
  }
  dz_[0] = dz_[1]; //hack
  zm_[0] = z_[0]-0.5*dz_[0];
  for(int j=0; j<num_points_-1; ++j) {
    dzm_[j] = zm_[j+1]-zm_[j];
  }
  dzm_[num_points_-1] = dz_[num_points_-1];

  // Copy global grid to local grid
  // boundary "ghost" nodes use boundary values
  for (int j=-nover_; j<num_local_points_+nover_; ++j) {
    int jglobal = j + my_pe_*num_local_points_;
    jglobal = std::max(0,jglobal);
    jglobal = std::min(jglobal, num_points_-1);
    dz_local_[nover_+j] = dz_[jglobal];
    dzm_local_[nover_+j] = dzm_[jglobal];
    inv_dz_local_[nover_+j] = 1.0/dz_[jglobal];
    inv_dzm_local_[nover_+j] = 1.0/dzm_[jglobal];
  }

  //set extended (with ghost cells) work arrays y_ext and rhs_ext
  y_ext_.assign( (num_local_points_+(2*nover_))*num_states_, 0.0);
  rhs_ext_.assign( (num_local_points_+(2*nover_))*num_states_, 0.0);

  rhsConv_.assign(num_local_points_*num_states_,0.0);

  // Create relative volume arrays
  rel_vol_.assign(num_local_points_,0.0);
  rel_vol_ext_.assign(num_local_points_+2*nover_,0.0);

}

void FlameParams::SetTfix() {
  // Set jfix and Tfix
  int local_jfix = num_points_+1;
  double local_Tfix = 0.0;
  double temperature_fix_target = 0.0;
  std::string tfix_setting = string_options_["temperature_fix_setting"];
  double tfix_value = double_options_["temperature_fix_value"];
  if(tfix_setting == "Delta") {
    temperature_fix_target = inlet_temperature_ + tfix_value/reference_temperature_;
  } else if(tfix_setting == "Absolute") {
    temperature_fix_target = tfix_value/reference_temperature_;
  } else if(tfix_setting == "InOutMix") {
    temperature_fix_target = tfix_value*inlet_temperature_ + (1-tfix_value)*outlet_temperature_;
  }
  for(int j=0; j<num_local_points_; ++j) {
    int jglobal = j + my_pe_*num_local_points_;
    int temp_id = j*num_states_+num_species_ + 1;
    if( y_[temp_id] > temperature_fix_target ) {
      local_jfix = jglobal;
      local_Tfix = y_[temp_id];
      break;
    }
  }
#ifdef ZERORK_MPI
  struct { int index; int rank; } local_fix_values, global_fix_values;
  local_fix_values.index = local_jfix;
  local_fix_values.rank = my_pe_;
  MPI_Allreduce(&local_fix_values ,&global_fix_values, 1, MPI_2INT, MPI_MINLOC, comm_);
  j_fix_ = global_fix_values.index;
  MPI_Bcast(&local_Tfix, 1, MPI_DOUBLE, global_fix_values.rank, comm_);
  temperature_fix_ = local_Tfix;
#else
  j_fix_ = local_jfix;
  temperature_fix_ = local_Tfix;
#endif
  if(j_fix_ == num_points_ + 1) {
    error_status_ = 6;
  }
}

// requires SetGrid to be set first
void FlameParams::SetMemory()
{
  const int num_nonzeros_zerod = reactor_->GetJacobianSize();
  const int num_reactors = num_local_points_;
  int num_local_states = num_local_points_*num_states_;
  int num_total_states = num_points_*num_states_;
  num_nonzeros_loc_ = 0;

  //std::vector<int> column_id;
  std::vector<int> row_id;
  int  dense_id;

  std::vector<int> row_id_zerod, col_id_zerod;

  mti_mf_.resize(num_species_,0.0); 
  mti_mf_grad_.resize(num_species_, 0.0);
  mti_temp_grad_.resize(1, 0.0); // num_dimensions
  mti_pres_grad_.resize(1, 0.0); // num_dimensions
  // Setup the MassTransportFactory data structure
  transport_input_.num_dimensions_        = 1;
  transport_input_.ld_grad_temperature_   = 1;
  transport_input_.ld_grad_pressure_      = 1;
  transport_input_.ld_grad_mass_fraction_ = num_species_;
  transport_input_.mass_fraction_         = &mti_mf_[0];
  transport_input_.grad_temperature_      = &mti_temp_grad_[0];
  transport_input_.grad_pressure_         = &mti_pres_grad_[0];
  transport_input_.grad_mass_fraction_    = &mti_mf_grad_[0];

  // Apply constant pressure approximation
  transport_input_.pressure_         = pressure_;
  transport_input_.grad_pressure_[0] = 0.0;

  // create and set the inverse molecular mass array
  inv_molecular_mass_.assign(num_species_, 0.0);
  reactor_->GetSpeciesMolecularWeight(&inv_molecular_mass_[0]);
  for(int j=0; j<num_species_; ++j) {
    inv_molecular_mass_[j] = 1.0/inv_molecular_mass_[j];
  }

  // "old" state vector for pseudo unsteady
  y_old_.assign(num_local_points_*num_states_, 0.0);

  // create the workspace for the species specific heats
  species_specific_heats_.assign(num_species_*(num_local_points_+1), 0.0);

  // create the workspace for the species mass fluxes and Lewis numbers
  species_mass_flux_.assign(num_species_*(num_local_points_+1), 0.0); //larger size for derivatives
  species_lewis_numbers_.assign((num_local_points_+1)*num_species_, 0.0);

  // create the workspace for the flux interface conductivities
  thermal_conductivity_.assign(num_local_points_+1, 0.0);//larger size for derivatives

  // create the workspace for the mixture specific heat at each grid point
  mixture_specific_heat_.assign(num_local_points_+1, 0.0);

  // create the workspace for the mid point mixture specific heat
  mixture_specific_heat_mid_.assign(num_local_points_+1, 0.0);

  // create the workspace for the mid point mixture molecular mass
  molecular_mass_mix_mid_.assign(num_local_points_+1, 0.0);

  // Get number of off-diagonals terms to keep in block jacobian
  num_off_diagonals_ = nover_*num_states_;

  // Use SuperLU serial
  superlu_serial_ = true;
  if(npes_ > 1) {
    superlu_serial_ = false;
  }

  // create the Jacobian data structures needed for the banded and SPGMR
  // integrators
  // Set Jacobian parameters
  // Case 2 is a block tridiagonal (potentially sparse) matrix solved with SuperLU
  if(integrator_type_ == 2) {
    row_id_zerod.assign(num_nonzeros_zerod, 0);
    col_id_zerod.assign(num_nonzeros_zerod, 0);

    // The Jacobian pattern is assumed to be in compressed column storage
    reactor_->GetJacobianPattern(&row_id_zerod[0],
                                 &col_id_zerod[0]);

    // N.B.: Currently integrator_type_ == 2 enforces dense blocks
    //       This may be reevaluated
#if 1
    // Set to 1 to force dense blocks
    dense_to_sparse_.assign(num_states_*num_states_, 1);
    dense_to_sparse_offdiag_.assign(num_states_*num_states_, 1);

#else
    dense_to_sparse_.assign(num_states_*num_states_, 0);
    dense_to_sparse_offdiag_.assign(num_states_*num_states_, 0);

    // Chemical jacobian pattern -- only for diagonal block
    for (int j=0; j<num_nonzeros_zerod; j++) {
      dense_id = num_states_*col_id_zerod[j] + row_id_zerod[j];
      dense_to_sparse_[dense_id] = 1;
      //dense_to_sparse_offdiag_[dense_id] = 1; //for off-diagonal too?
    }

    for (int j=0; j<num_states_; j++) {
      for (int i=0; i<num_states_; i++) {
	dense_id = num_states_*j + i;

	//Dense rows and columns for local mdot
	if(j==num_states_-2 || i==num_states_-2) {
	  //dense_to_sparse_[dense_id] = 1; //already in getjacobianpattern
	}
	//Dense rows and columns for local T
	if(j==num_states_-1 || i==num_states_-1) {
	  //dense_to_sparse_[dense_id] = 1; //already in getjacobianpattern
	}
	//Dense rows for off-diagonal T?
	if(i==num_states_-1) {
	  //dense_to_sparse_offdiag_[dense_id] = 1;
	}

	//Diagonal
	if(i==j) {
	  dense_to_sparse_[dense_id] = 1;
	  dense_to_sparse_offdiag_[dense_id] = 1;
	}
      }
    }
#endif

    // Count non-zeros
    int width = 2*num_off_diagonals_ + 1;
    int i1, i2;
    for (int group=1; group<=width; group++) {
      for (int j=group-1; j<num_local_points_*num_states_; j+=width) {
	int jglobal = j + my_pe_*num_local_points_*num_states_;
	int jstate = jglobal % num_states_;
	i1 = std::max(0, jglobal - jstate - num_states_);
	i2 = std::min(jglobal + (num_states_-1 - jstate) + num_states_, num_points_*num_states_-1);
	for (int i=i1; i<=i2; i++) {
	  int istate = i % num_states_;
	  dense_id = num_states_*jstate + istate; //here j is column and i is row
	  //Diagonal block
	  if (i>= jglobal-jstate && i<=jglobal+num_states_-1-jstate) {
	    if (dense_to_sparse_[dense_id] == 1)
	      num_nonzeros_loc_++;
	  }
	  //Off-diagonal block
	  if (i<jglobal-jstate || i>jglobal+num_states_-1-jstate) {
	    if (dense_to_sparse_offdiag_[dense_id] == 1)
	      num_nonzeros_loc_++;
	  }
	}
      }
    }

    // SuperLU
    reactor_jacobian_dist_.assign(num_nonzeros_loc_, 0.0);
    col_id_.assign(num_nonzeros_loc_, 0);
    row_id.assign(num_nonzeros_loc_, 0);
    row_sum_.assign(num_local_points_*num_states_+1,0);

    // Get pattern "manually" for now
    int innz=0;
    // For compressed row storage j is the row and i is the column
    for (int j=0; j<num_local_points_*num_states_; j++) {
      int jglobal = j + my_pe_*num_local_points_*num_states_;
      int jstate = jglobal % num_states_;
      i1 = std::max(0, jglobal - jstate - num_states_);
      i2 = std::min(jglobal + (num_states_-1 - jstate) + num_states_, num_points_*num_states_-1);
      for (int i=i1; i<=i2; i++) {
	// Compressed column storage
	// row_id_[innz] = i;
	// column_id[innz] = j;

	// Compressed row storage
	int istate = i % num_states_;
	dense_id = num_states_*istate + jstate; //i is column and j is row

	//Diagonal block. Should use GetJacobianPattern?
	if (i>= jglobal-jstate && i<=jglobal+num_states_-1-jstate) {
	  if (dense_to_sparse_[dense_id] == 1) {
	    row_id[innz] = j;
	    col_id_[innz] = i;
	    innz++;
	  }
	}
	//Off-diagonal block.
	if (i<jglobal-jstate || i>jglobal+num_states_-1-jstate) {
	  if (dense_to_sparse_offdiag_[dense_id] == 1) {
	    row_id[innz] = j;
	    col_id_[innz] = i;
	    innz++;
	  }
	}
      } // for i=i1 to i2
    } //for j=0 to < num_local_states

    for(int j=0; j<num_nonzeros_loc_; ++j) {
        ++row_sum_[row_id[j]+1];
    } // for(int j=0; j<num_nonzeros; ++j)

    // At this point the column sum contains just the count per column, need
    // to loop through to change the count to a sum
    for(int j=0; j<num_local_points_*num_states_; ++j) {
      row_sum_[j+1] += row_sum_[j];
    }

    if(superlu_serial_) {
      sparse_matrix_ = new SparseMatrix(num_local_points_*num_states_, num_nonzeros_loc_);
      if(sparse_matrix_ == NULL) {
	error_status_ = 3;
	return;
      }
#ifdef ZERORK_MPI
    } else {
      sparse_matrix_dist_ = new SparseMatrix_dist(num_local_points_*num_states_, num_nonzeros_loc_, comm_);
      if(sparse_matrix_dist_ == NULL) {
	error_status_ = 3;
	return;
      }
#endif
    }

    if(store_jacobian_) {
      saved_jacobian_dist_.assign(num_nonzeros_loc_, 0.0);
    }
  } // if integrator_type == 2

  // Case 3 is an approximately factorized matrix
  // local sparse matrices at each grid point for the chemical Jacobian (SuperLU)
  // global tridiagonal matrix for the transport Jacobian (LAPACK)
  if(integrator_type_ == 3) {
    //SuperLU
    row_id_chem_.assign(num_nonzeros_zerod, 0);
    column_id_chem_.assign(num_nonzeros_zerod, 0);
    column_sum_chem_.assign(num_states_+1,0);

    reactor_jacobian_chem_.assign(num_nonzeros_zerod, 0.0);

    // The Jacobian pattern is assumed to be in compressed column storage
    reactor_->GetJacobianPattern(&row_id_chem_[0],
				 &column_id_chem_[0]);

    for(int j=0; j<num_nonzeros_zerod; ++j) {
      ++column_sum_chem_[column_id_chem_[j]+1];
    }

    for(int j=0; j<num_states_; ++j) {
      column_sum_chem_[j+1] += column_sum_chem_[j];
    }


    diagonal_id_chem_.assign(num_states_, -1); // use negative one as flag to
    // check if there are any missing
    // diagonal terms that would prevent
    // its use for the sparse factorization
    for(int j=0; j<num_nonzeros_zerod; ++j) {
      if(row_id_chem_[j] == column_id_chem_[j]) {
	diagonal_id_chem_[row_id_chem_[j]] = j;
      }
    }
    for(int j=0; j<num_states_; ++j) {
      if(diagonal_id_chem_[j] == -1) {
	valid_jacobian_structure_ = false;
        error_status_ = 4;
        return;
      }
    }

    sparse_matrix_chem_.assign(num_reactors, NULL);
    for(int j=0; j<num_reactors; ++j) {
      sparse_matrix_chem_[j] = new SparseMatrix(num_states_, num_nonzeros_zerod);
      if(sparse_matrix_chem_[j] == NULL) {
        error_status_ = 4;
        return;
      }
    }

    if(store_jacobian_) {
      saved_jacobian_chem_.assign(num_nonzeros_zerod*num_reactors, 0.0);
    }

    // LAPACK
    int storage = 4*1 + 1;
    banded_jacobian_.assign(num_local_points_*storage * num_states_, 0.0);

    // test over whole domain
    num_states_per_proc_ = (num_states_ + npes_ - 1)/npes_;
    if(num_states_ - (npes_-1)*num_states_per_proc_ < 1) {
      // "Too few states per processor. Try using fewer processors.\n";
      error_status_ = 5;
      return;
    }
    if(my_pe_ == npes_-1) {
      num_states_local_ = num_states_ - (npes_-1)*num_states_per_proc_;
    } else {
      num_states_local_ = num_states_per_proc_;
    }
    banded_jacobian2_.assign(num_points_*storage * num_states_local_, 0.0);
    banded_jacobian_serial_.assign(num_points_*4*num_states_local_, 0.0);
    pivots_serial_.assign(num_points_*num_states_local_, 0);

  } // if integrator_type == 3
}

void FlameParams::GetTemperatureAndMassFractions(double* T, double* mass_fractions)
{
#ifdef ZERORK_MPI
  //Collect output
  int remainder = num_points_ % npes_;
  if(remainder != 0) {
    error_status_ = 1;
    return;
  }
  if(my_pe_ < remainder) {
    num_local_points_ += 1;
  }
  std::vector<int> gather_counts(npes_,0);
  std::vector<int> gather_counts_mf(npes_,0);
  std::vector<int> gather_displs(npes_,0);
  std::vector<int> gather_displs_mf(npes_,0);
  for(int j = 0; j < npes_; ++j) {
    gather_counts[j] = num_points_/npes_;
    if(j<remainder) gather_counts[j] += 1;
    gather_counts_mf[j] = gather_counts[j]*num_species_;
    if(j>0) {
      gather_displs[j] = gather_displs[j-1]+gather_counts[j-1];
      gather_displs_mf[j] = gather_displs_mf[j-1]+gather_counts_mf[j-1];
    }
  }

  std::vector<double> mass_fractions_local(num_local_points_*num_species_,0.0);
  std::vector<double> T_local(num_local_points_,0.0);
  for(int j = 0; j < num_local_points_; ++j) {
    for(int k = 0; k < num_species_; ++k) {
      mass_fractions_local[j*num_species_ + k] = y_[j*num_states_ +k];
    }
    T_local[j] = y_[(j+1)*num_states_ - 1]*reference_temperature_;
  }
  if(my_pe_ == 0) {
    for(int k = 0; k < num_species_; ++k) {
      mass_fractions[k] = inlet_mass_fractions_[k];
    }
    T[0] = inlet_temperature_*reference_temperature_;
  }     

  MPI_Gatherv(mass_fractions_local.data(), gather_counts_mf[my_pe_], MPI_DOUBLE,
              &mass_fractions[num_species_], gather_counts_mf.data(), gather_displs_mf.data(), MPI_DOUBLE, 0, comm_);
  MPI_Gatherv(T_local.data(), gather_counts[my_pe_], MPI_DOUBLE,
              &T[1], gather_counts.data(), gather_displs.data(), MPI_DOUBLE, 0, comm_);
#else
  for(int j = 0; j < num_local_points_; ++j) {
    for(int k = 0; k < num_species_; ++k) {
      mass_fractions[(j+1)*num_species_ + k] = y_[j*num_states_ +k];
    }
    T[j+1] = y_[(j+1)*num_states_ - 1]*reference_temperature_;
  }
  for(int k = 0; k < num_species_; ++k) {
    mass_fractions[k] = inlet_mass_fractions_[k];
  }
  T[0] = inlet_temperature_*reference_temperature_;
#endif
} // void FlameParams::GetTemperatureAndMassFractions()

