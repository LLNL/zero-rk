#include <stdlib.h> // exit()

#include <string>   // C++ std::string

#include "reactor_base.h"
#include "variable_volume_reactor.h"


// implementation class for the VariableVolumeReactor
class VariableVolumeReactor::Impl: public ReactorBase
{
 public:
  // member functions that must be defined for the ReactorBase class
  Impl(const char mechanism_name[],
       const char thermodynamics_name[],
       const char parser_log_name[],
       MatrixType matrix_type,
       Volume *volume_function,
       HeatLoss *heat_loss_function);
  ~Impl();
  ReactorError GetTimeDerivative(const double reactor_time,
                                 const double state[],
				 double derivative[]);

  ReactorError GetChemicalHeatReleaseRate(const double reactor_time,
                                          const double state[],
				          double *heat_release_rate);

  ReactorError GetJacobian(const double reactor_time,
                           const double state[],
			   double jacobian[]);

  int GetJacobianPattern(int row_id[],
                         int column_id[]);
  // --------------------------------------
  // additional functions specific to the reactor type
  void SetReferenceScales(const double ref_moles,
                          const double ref_temperature);
  void GetReferenceScales(double *ref_moles,
                          double *ref_temperature) const;

  void SetVolumeMultiplier(const double v_multiplier) 
  {v_multiplier_ = v_multiplier;}
  double GetVolumeMultiplier() const {return v_multiplier_;}

  // specific heat calculation
  double GetMixtureSpecificHeat_Cv(const double state[]);

 private:
  int BuildSparseJacobianArrays();
  ReactorError GetSparseJacobian(const double reactor_time,
                                 const double state[],
			         double jacobian[]);
  void GetJacobianColumnFromPerturbation(const int column_id,
                                         const double dstate,
                                         const double reactor_time,
                                         const double original_state[],
                                         const double original_derivative[],
                                         double perturbed_state[],
                                         double perturbed_derivative[],
                                         double jacobian_column[]);

  Volume *volume_function_;
  HeatLoss *heat_loss_function_;
  double ref_moles_;
  double inv_ref_moles_;
  double ref_temperature_;
  double inv_ref_temperature_;
  bool use_scaled_state_;

  double v_multiplier_;

  // time derivative storage arrays
  std::vector<double> concentrations_;
  std::vector<double> internal_energies_;
  std::vector<double> specific_heats_;
  std::vector<double> net_reaction_rates_;
  std::vector<double> creation_rates_;
  std::vector<double> destruction_rates_;
  std::vector<double> step_rates_;

  // Jacobian storage arrays
  std::vector<int> destroy_concentration_id_;
  std::vector<int> destroy_step_id_;
  std::vector<int> destroy_sparse_id_;
  std::vector<int> create_concentration_id_;
  std::vector<int> create_step_id_;
  std::vector<int> create_sparse_id_;
  std::vector<int> jacobian_row_id_;
  std::vector<int> jacobian_column_id_;
  std::vector<int> jacobian_column_sum_;
  // Jacobian calculation work arrays
  std::vector<double> inv_concentrations_;
  std::vector<double> original_state_;
  std::vector<double> perturbed_state_;
  std::vector<double> original_derivative_;
  std::vector<double> perturbed_derivative_;

};

// Implementation of the pure member functions from the ReactorBase class
VariableVolumeReactor::Impl::Impl(const char mechanism_name[],
                                  const char thermodynamics_name[],
                                  const char parser_log_name[],
                                  MatrixType matrix_type,
                                  Volume *volume_function,
                                  HeatLoss *heat_loss_function)
{
  int jacobian_size;
  std::string info;
  std::string current_name;
  std::string species_prefix;
  std::vector<std::string> state_names;
  zerork::mechanism *mechanism_ptr;


  // set constructor arguments
  SetMatrixType(matrix_type);
  volume_function_ = volume_function;
  heat_loss_function_ = heat_loss_function;

  // set reference scales for state quanitites
  use_scaled_state_ = false;
  inv_ref_moles_ = ref_moles_ = 1.0;
  inv_ref_temperature_ = ref_temperature_ = 1.0;
  // set the initial volume multiplier
  v_multiplier_ = 1.0;
 
  BuildMechanism(mechanism_name,
                 thermodynamics_name,
                 parser_log_name);

  // create the vector of state names
  state_names.clear();
  mechanism_ptr = GetMechanism();
  species_prefix = std::string("Moles_");
  const int num_species = mechanism_ptr->getNumSpecies();
  const int num_steps   = mechanism_ptr->getNumSteps();

  for(int j=0; j<num_species; ++j) {
    current_name = std::string(mechanism_ptr->getSpeciesName(j));
    current_name.insert(0,species_prefix); // prepend prefix
    state_names.push_back(current_name);
  }
  state_names.push_back(std::string("TotalMoles"));
  state_names.push_back(std::string("Temperature"));
  BuildStateNamesMap(state_names);
  const int num_states = static_cast<int>(state_names.size());
  //printf("state_name.size() = %lu\n",state_names.size());
  //fflush(stdout);
  //printf("BuildStateNamesMap(state_names) = %d\n",
  //       BuildStateNamesMap(state_names));
  //printf("GetNumStates() = %d\n",GetNumStates());
  //fflush(stdout);

  // set the jacobian size with the number of non-zero terms
  // and build the associated data structures for the sparse Jacobian
  // if used.
  if(matrix_type == DENSE_COL_MAJOR) {
    jacobian_size  = GetNumStates();
    jacobian_size *= jacobian_size;
  } else if(matrix_type == COMPRESSED_COL_STORAGE) {
    jacobian_size = BuildSparseJacobianArrays();
  } else {
    printf("ERROR: unrecognized matrix_type %d\n",
           static_cast<int>(matrix_type));
    fflush(stdout);
    exit(-1);
  }
  SetJacobianSize(jacobian_size);

  // provide basic info about the reactor model
  info  = "Variable volume reactor\n";
  info += "with number of species N = ";
  // TO DO: make the string formation cleaner
  char tmp[32];
  sprintf(tmp,"%d",num_species);
  info += std::string(tmp);
  info += std::string("\nstate[0]:   ") + state_names[0];
  info += std::string("\nstate[N-1]: ") + state_names[num_species-1];
  info += std::string("\nstate[N]:   ") + state_names[num_species];
  info += std::string("\nstate[N+1]: ") + state_names[num_species+1];
  info += "\n-----------------------\n";
  SetReactorInfo(info); 

  // pre-assign the size of the internal time-derivative vectors
  concentrations_.assign(num_species,0.0);
  internal_energies_.assign(num_species,0.0);
  specific_heats_.assign(num_species,0.0);
  net_reaction_rates_.assign(num_species,0.0);
  creation_rates_.assign(num_species,0.0);
  destruction_rates_.assign(num_species,0.0);
  step_rates_.assign(num_steps,0.0);
  // pre-assign the size of the internal jacobian space vectors
  inv_concentrations_.assign(num_species,0.0);
  original_state_.assign(num_states,0.0);
  perturbed_state_.assign(num_states,0.0);
  original_derivative_.assign(num_states,0.0);
  perturbed_derivative_.assign(num_states,0.0);
}

VariableVolumeReactor::Impl::~Impl() 
{
  DestroyMechanism();
}

ReactorError
  VariableVolumeReactor::Impl::GetTimeDerivative(const double reactor_time,
                                                 const double state[],
                                                 double derivative[])
{
  // define local constants of class members to enable loop vectorization
  const int num_species            = GetNumSpecies();
  const double ref_moles           = ref_moles_;
  const double ref_temperature     = ref_temperature_;
  const double inv_ref_moles       = inv_ref_moles_;
  const double inv_ref_temperature = inv_ref_temperature_;
  const double total_moles         = ref_moles_*state[num_species];
  const double temperature         = ref_temperature*state[num_species+1];
  
  double volume, dvolume_dt, heat_loss, molar_cv_mix, RuT;
  double mole_sum, internal_energy_sum;
  zerork::mechanism *mechanism_ptr;

  mechanism_ptr = GetMechanism();

  // compute the product of the universal gas constant and temperature
  RuT = temperature*mechanism_ptr->getGasConstant();

  volume_function_->GetVolume(reactor_time,
                              state,
                              &volume,
                              &dvolume_dt);
  volume     *= v_multiplier_;
  dvolume_dt *= v_multiplier_; 

  heat_loss_function_->GetHeatLoss(reactor_time,
                                   state,
				   &heat_loss);
  // compute concentration, total_moles and temperature from the current state
  for(int j=0; j<num_species; ++j) {
    concentrations_[j] = ref_moles*state[j]/volume;
  }
  // compute the rate of change of the species concentration
  mechanism_ptr->getReactionRates_perturbROP(temperature,
                                             &concentrations_.at(0),
                                             GetAMultipliers(),
                                             &net_reaction_rates_.at(0),
                                             &creation_rates_.at(0),
                                             &destruction_rates_.at(0),
                                             &step_rates_.at(0));

  // compute the molar specific heat of the mixture
  molar_cv_mix = mechanism_ptr->getMolarCvFromTC(temperature,
                                                 &concentrations_.at(0),
                                                 &specific_heats_.at(0));
  // compute the species internal energies
  mechanism_ptr->getIntEnergy_RT(temperature,
                                 &internal_energies_.at(0));

  // compute the rate of change of the moles of each species
  mole_sum = 0.0;
  internal_energy_sum = 0.0;
  for(int j=0; j<num_species; ++j) {
    derivative[j] = net_reaction_rates_[j]*volume;
    mole_sum += derivative[j];
    internal_energy_sum += derivative[j]*internal_energies_[j];
  }

  // compute the rate of temperature change
  // dT/dt = -1/(Cv_mix*n_t)
  derivative[num_species+1] = -(1.0/(molar_cv_mix*total_moles))*
    (internal_energy_sum*RuT + heat_loss + 
     dvolume_dt*total_moles*RuT/volume);

  // rescale the state vector derivative
  for(int j=0; j<num_species; ++j) {
    derivative[j] *= inv_ref_moles_;
  }
  // rate of change of the total moles
  derivative[num_species]    = inv_ref_moles*mole_sum; 
  derivative[num_species+1] *= inv_ref_temperature;
   
  return NONE; // no error
}
ReactorError
  VariableVolumeReactor::Impl::GetChemicalHeatReleaseRate(const double reactor_time, const double state[], double *heat_release_rate)
{
  // define local constants of class members to enable loop vectorization
  const int num_species            = GetNumSpecies();
  const double ref_moles           = ref_moles_;
  const double ref_temperature     = ref_temperature_;
  const double temperature         = ref_temperature*state[num_species+1];
  
  double RuT, volume, dvolume_dt;
  double internal_energy_sum;
  zerork::mechanism *mechanism_ptr;

  mechanism_ptr = GetMechanism();

  // compute the product of the universal gas constant and temperature
  RuT = temperature*mechanism_ptr->getGasConstant();


  volume_function_->GetVolume(reactor_time,
                              state,
                              &volume,
                              &dvolume_dt);
  volume     *= v_multiplier_;
  dvolume_dt *= v_multiplier_; 

  // compute concentration, total_moles and temperature from the current state
  for(int j=0; j<num_species; ++j) {
    concentrations_[j] = ref_moles*state[j]/volume;
  }
  // compute the rate of change of the species concentration
  mechanism_ptr->getReactionRates_perturbROP(temperature,
                                             &concentrations_.at(0),
                                             GetAMultipliers(),
                                             &net_reaction_rates_.at(0),
                                             &creation_rates_.at(0),
                                             &destruction_rates_.at(0),
                                             &step_rates_.at(0));

  // compute the species internal energies
  mechanism_ptr->getIntEnergy_RT(temperature,
                                 &internal_energies_.at(0));

  // compute the rate of change of the moles of each species
  internal_energy_sum = 0.0;
  for(int j=0; j<num_species; ++j) {
    internal_energy_sum += net_reaction_rates_[j]*volume*internal_energies_[j];
  }
  // Change sign of chemical heat release rate to be consistent with Heywood's
  // definition in eqs. 9.24 and 9.25
  (*heat_release_rate) = -internal_energy_sum*RuT;  
   
  return NONE; // no error
}

ReactorError
  VariableVolumeReactor::Impl::GetJacobian(const double reactor_time,
                                           const double state[],
                                           double jacobian[])
{

  if(GetMatrixType() == DENSE_COL_MAJOR) {
    printf("ERROR: GetJacobian not currently defined for\n");
    printf("       DENSE_COL_MAJOR.\n");
    return UNKNOWN;
  } else if(GetMatrixType() == COMPRESSED_COL_STORAGE) {
    return GetSparseJacobian(reactor_time, state, jacobian);
  }

  return UNKNOWN;
}
int VariableVolumeReactor::Impl::GetJacobianPattern(int row_id[],
                                                    int column_id[])
{
  const int num_states    = GetNumStates();
  const int jacobian_size = GetJacobianSize(); 
  if(row_id != NULL && column_id != NULL) {

    if(GetMatrixType() == DENSE_COL_MAJOR) {
      int dense_id = 0;
      for(int j=0; j<num_states; ++j) {
        // column - j
        for(int k=0; k<num_states; ++k) {
          // row - k
          row_id[dense_id]    = k;
          column_id[dense_id] = j;
          ++dense_id;
        }
      }

    } else if(GetMatrixType() == COMPRESSED_COL_STORAGE) {
      for(int j=0; j<jacobian_size; ++j) {
        row_id[j]    = jacobian_row_id_[j];
        column_id[j] = jacobian_column_id_[j];
      }
    }
  }

  return jacobian_size;
}

int VariableVolumeReactor::Impl::BuildSparseJacobianArrays()
{
  std::map<int, int> dense_to_sparse_map;
  std::map<int, int>::iterator map_iter;
  zerork::mechanism *mechanism_ptr;
  int row_id, column_id, dense_id, sparse_id;


  mechanism_ptr = GetMechanism();
  if(mechanism_ptr == NULL) {
    printf("ERROR: In BuildSparseJacobianArrays(),\n");
    printf("       mechanism is not allocated.\n");
    fflush(stdout);
    return -1;
  }
  const int num_species = mechanism_ptr->getNumSpecies();
  const int num_steps   = mechanism_ptr->getNumSteps();
  const int num_states  = num_species + 2; // total moles & temperature

  // clear the jacobian arrays
  destroy_concentration_id_.clear();
  destroy_step_id_.clear();
  destroy_sparse_id_.clear();
  create_concentration_id_.clear();
  create_step_id_.clear();
  create_sparse_id_.clear();
  jacobian_row_id_.clear();
  jacobian_column_id_.clear();
  jacobian_column_sum_.clear();
  // temporary map to store distinct array points
  // also will count the total number of jacobian terms (elementary steps)
  // that w
  dense_to_sparse_map.clear();
  // 
  for(int j=0; j<num_steps; ++j) {

    int num_reactants = mechanism_ptr->getOrderOfStep(j);
    int num_products  = mechanism_ptr->getNumProductsOfStep(j);

    for(int k=0; k<num_reactants; ++k) {
      // get the index of the species being perturbed for the Jacobian
      column_id = mechanism_ptr->getSpecIdxOfStepReactant(j,k);
      
      // set the Jacobian terms related to the destruction of a particular
      // species affected by the perturbation of species column_id
      for(int m=0; m<num_reactants; ++m) {

        row_id   = mechanism_ptr->getSpecIdxOfStepReactant(j,m);
        dense_id = row_id+column_id*num_states;

        destroy_concentration_id_.push_back(column_id);
        destroy_step_id_.push_back(j);
        destroy_sparse_id_.push_back(dense_id);

        // record position in dense matrix
        map_iter = dense_to_sparse_map.find(dense_id);
        if(map_iter == dense_to_sparse_map.end()) {
          // first term for dense_id
          dense_to_sparse_map[dense_id] = 1;
        } else {
          // increment number of terms for dense_id
          dense_to_sparse_map[dense_id] += 1; 
        }

      } // m-loop over reactants of step j

      // set the Jacobian terms related to the creation of a particular
      // species affected by the perturbation of species column_id
      for(int m=0; m<num_products; ++m) {
        row_id = mechanism_ptr->getSpecIdxOfStepProduct(j,m);
        dense_id = row_id+column_id*num_states;

        create_concentration_id_.push_back(column_id);
        create_step_id_.push_back(j);
        create_sparse_id_.push_back(dense_id);

        // record position in dense matrix
        map_iter = dense_to_sparse_map.find(dense_id);
        if(map_iter == dense_to_sparse_map.end()) {
          // first term for dense_id
          dense_to_sparse_map[dense_id] = 1;
        } else {
          // increment number of terms for dense_id
          dense_to_sparse_map[dense_id] += 1; 
        }

      } // m-loop over products of step j

    } // k-loop over reactants of step j

  } // j-loop over steps
  //printf("# Number of distinct elements in composition jacobian: %d\n",
  //       static_cast<int>(dense_to_sparse_map.size()));
  //printf("# Number of jacobian terms for elementary steps      : %d (create)\n",
  //        static_cast<int>(create_sparse_id_.size()));
  //printf("# Number of jacobian terms for elementary steps      : %d (destroy)\n",
  //       static_cast<int>(destroy_sparse_id_.size()));

  // add dense rows for the total moles and temperature states
  // affected by perturbations in the species concentrations
  for(int j=0; j<num_species; ++j) {
    // add total moles row (row_id = num_species)
    dense_id = num_species + j*num_states;
    dense_to_sparse_map[dense_id] = 1;
    // add temperature row (row_id = num_species+1)
    dense_id = num_species + 1 + j*num_states;
    dense_to_sparse_map[dense_id] = 1;
  }
  // add dense columns for all the states affected by perturbations in
  // the total moles and the temperature
  for(int j=0; j<num_states; ++j) {
    // add total moles column (column_id = num_species)
    dense_id = j + num_species*num_states;
    dense_to_sparse_map[dense_id] = 1;
    // add temperature column (column_id = num_species+1)
    dense_id = j + (num_species+1)*num_states;
    dense_to_sparse_map[dense_id] = 1;
  }

  // add the diagonal for all states
  for(int j=0; j<num_states; ++j) {
    dense_id = j+j*num_states;
    // record position in dense matrix
    map_iter = dense_to_sparse_map.find(dense_id);
    if(map_iter == dense_to_sparse_map.end()) {
      // first term for dense_id
      dense_to_sparse_map[dense_id] = 1;
    }
  }

  // record the non-zero addresses for the jacobian_column_id_ and
  // jacobian_row_id_, and update the dense_to_sparse_map storing the sparse_id
  // for all the existing non-zero dense elements
  dense_id = 0;
  sparse_id = 0;
  for(int j=0; j<num_states; ++j) {
    // column-j
    jacobian_column_sum_.push_back(sparse_id);
    for(int k=0; k<num_states; ++k) {
      // row-k
      map_iter = dense_to_sparse_map.find(dense_id);
      if(map_iter != dense_to_sparse_map.end()) {
        // found non-zero element recorded in dense_to_sparse_map
        dense_to_sparse_map[dense_id] = sparse_id;
        jacobian_row_id_.push_back(k);
        jacobian_column_id_.push_back(j);
        ++sparse_id;
      }
      ++dense_id;
    }  
  }
  jacobian_column_sum_.push_back(sparse_id);

  // use the dense to sparse map to transform the dense_id originally stored
  // in destroy_sparse_id_ and create_sparse_id_.
  const int num_destroy = static_cast<int>(destroy_sparse_id_.size());
  for(int j=0; j<num_destroy; ++j) {
    dense_id = destroy_sparse_id_[j];
    destroy_sparse_id_[j] = dense_to_sparse_map[dense_id];
  }
  const int num_create  = static_cast<int>(create_sparse_id_.size());
  for(int j=0; j<num_create; ++j) {
    dense_id = create_sparse_id_[j];
    create_sparse_id_[j] = dense_to_sparse_map[dense_id];
  }

  return static_cast<int>(dense_to_sparse_map.size()); // number of non-zeros
}


ReactorError 
  VariableVolumeReactor::Impl::GetSparseJacobian(const double reactor_time,
                                                 const double state[],
			                         double jacobian[])
{
  const double perturb_factor = 1.0e-8;
  // define local constants of class members to enable loop vectorization
  const int num_species            = GetNumSpecies();
  //const int num_states             = GetNumStates();
  const int num_nonzeros           = GetJacobianSize();
  const double ref_moles           = ref_moles_;
  const double ref_temperature     = ref_temperature_;
  const double inv_ref_moles       = inv_ref_moles_;
  const double inv_ref_temperature = inv_ref_temperature_;
  //const double total_moles         = ref_moles_*state[num_species];
  const double temperature         = ref_temperature*state[num_species+1];
  const int num_destroy  = static_cast<int>(destroy_sparse_id_.size());
  const int num_create   = static_cast<int>(create_sparse_id_.size());
  
  double volume, dvolume_dt, heat_loss, total_cv, internal_energy_sum;
  double d_temperature, d_total_moles;
  double min_concentration;
  zerork::mechanism *mechanism_ptr;

  mechanism_ptr = GetMechanism();

  // initialize the jacobian array to zero
  for(int j=0; j<num_nonzeros; ++j) {
    jacobian[j] = 0.0;
  }


  volume_function_->GetVolume(reactor_time,
                              state,
                              &volume,
                              &dvolume_dt);
  volume     *= v_multiplier_;
  dvolume_dt *= v_multiplier_; 

  heat_loss_function_->GetHeatLoss(reactor_time,
                                   state,
				   &heat_loss);

  // one molecule in the volume is taken to be the minimum concentration
  original_state_[num_species]   = 0.0; // recompute the total moles
  original_state_[num_species+1] = state[num_species+1];
  min_concentration = 1.0/(mechanism_ptr->getAvogadroNumber()*volume);
  for(int j=0; j<num_species; ++j) {
    concentrations_[j] = ref_moles*state[j]/volume;
    if(concentrations_[j] < min_concentration) {
      concentrations_[j] = min_concentration;
    }
    original_state_[j] = concentrations_[j]*volume*inv_ref_moles;
    original_state_[num_species] += original_state_[j];
    inv_concentrations_[j] = 1.0/concentrations_[j];
  }

  // compute the rate of change of the species concentration
  mechanism_ptr->getReactionRates_perturbROP(temperature,
                                             &concentrations_[0],
                                             GetAMultipliers(),
                                             &net_reaction_rates_[0],
                                             &creation_rates_[0],
                                             &destruction_rates_[0],
                                             &step_rates_[0]);

  // use the step rates and inverse concentrations with the elementary
  // Jacobian term lists to compute dwdot[i]/dC[j]
  // Destruction terms
  for(int j=0; j<num_destroy; ++j) {
    int conc_id   = destroy_concentration_id_[j];
    int step_id   = destroy_step_id_[j];
    int sparse_id = destroy_sparse_id_[j];

    jacobian[sparse_id] -= step_rates_[step_id]*inv_concentrations_[conc_id];
  }
  // Creation terms
  for(int j=0; j<num_create; ++j) {
    int conc_id   = create_concentration_id_[j];
    int step_id   = create_step_id_[j];
    int sparse_id = create_sparse_id_[j];

    jacobian[sparse_id] += step_rates_[step_id]*inv_concentrations_[conc_id];
  }

  mechanism_ptr->getCv_R_IntEnergy_RT(temperature,
                                      &specific_heats_[0],
                                      &internal_energies_[0]);
  
  total_cv = 0.0;
  internal_energy_sum = 0.0;
  for(int j=0; j<num_species; ++j) {
    total_cv += concentrations_[j]*specific_heats_[j];
    internal_energy_sum += internal_energies_[j]*net_reaction_rates_[j];
  }
  total_cv *= volume; // sum (species moles * Cv/Ru)
  // sum (species molar production rate * species internal energy per mole)
  internal_energy_sum *= volume*temperature;

  // build the derivative correpsonding to the original state
  original_derivative_[num_species] = 0.0;
  for(int j=0; j<num_species; ++j) {
    original_derivative_[j] = net_reaction_rates_[j] * volume * inv_ref_moles;
    original_derivative_[num_species] += original_derivative_[j];
  }
  // unscaled temperature derivative
  original_derivative_[num_species+1] = -(1.0/total_cv)*
    (internal_energy_sum + heat_loss/mechanism_ptr->getGasConstant() +
     (dvolume_dt/volume)*temperature*original_state_[num_species]*ref_moles);
  // scaled temperature derivative
  original_derivative_[num_species+1] *= inv_ref_temperature;

  // reuse the internal energy sum to store the multiplicative factor for 
  // each species Cv/Ru
  internal_energy_sum = (internal_energy_sum+heat_loss)/(total_cv*total_cv);
  // convert multiplicative factor from [K/kmol/s] to [1/s]
  internal_energy_sum *= inv_ref_temperature*ref_moles;


  // compute the original 
  // TODO: make the Jacobian construction more flexible
  // assumes compressed column storage with dense rows and columns for
  // total moles and temperature
  for(int j=0; j<num_species; ++j) {
    int num_nonzero_rows = jacobian_column_sum_[j+1] - jacobian_column_sum_[j];
    int sparse_id = jacobian_column_sum_[j];
    double jacobian_sum = 0.0;
    double jacobian_energy_sum = 0.0;
    for(int k=0; k<(num_nonzero_rows-2); ++k) {
      // sum use for the total moles row
      jacobian_sum += jacobian[sparse_id];

      // sum used for the temperature row
      int row_id = jacobian_row_id_[sparse_id];
      jacobian_energy_sum += jacobian[sparse_id]*internal_energies_[row_id];

      ++sparse_id;
    }
    jacobian[sparse_id]   = jacobian_sum;
    jacobian[sparse_id+1] = -inv_ref_temperature*ref_moles*temperature*
      jacobian_energy_sum/total_cv + internal_energy_sum*specific_heats_[j];
    
  }

  // perturb the total moles
  d_total_moles = original_state_[num_species]*perturb_factor;
  GetJacobianColumnFromPerturbation(num_species,
                                 d_total_moles,
                                 reactor_time,
                                 &original_state_[0],
                                 &original_derivative_[0],
                                 &perturbed_state_[0],
                                 &perturbed_derivative_[0],
				 &jacobian[jacobian_column_sum_[num_species]]);
  // perturb the temperature
  d_temperature = original_state_[num_species+1]*perturb_factor;
  GetJacobianColumnFromPerturbation(num_species+1,
                               d_temperature,
                               reactor_time,
                               &original_state_[0],
                               &original_derivative_[0],
                               &perturbed_state_[0],
                               &perturbed_derivative_[0],
			       &jacobian[jacobian_column_sum_[num_species+1]]);


 return NONE;
}


void VariableVolumeReactor::Impl::
  GetJacobianColumnFromPerturbation(const int column_id,
                                    const double dstate,
                                    const double reactor_time,
                                    const double original_state[],
                                    const double original_derivative[],
                                    double perturbed_state[],
                                    double perturbed_derivative[],
                                    double jacobian_column[])
{
  const int num_states = GetNumStates();
  double dstate_exact;
  // copy the orginal state
  for(int j=0; j<num_states; ++j) {
    perturbed_state[j] = original_state[j];
  }
  perturbed_state[column_id] += dstate;
  // compute the perturbation as represented in the floating point arithmetic
  dstate_exact = perturbed_state[column_id]-original_state[column_id];

  // compute the time derivative at the perturbed state
  GetTimeDerivative(reactor_time,perturbed_state,perturbed_derivative);

  // compute the jacobian column using a finite difference approximation
  dstate_exact = 1.0/dstate_exact;
  for(int j=0; j<num_states; ++j) {
    jacobian_column[j] =
      dstate_exact*(perturbed_derivative[j]-original_derivative[j]);
  }
}



void VariableVolumeReactor::Impl::SetReferenceScales(const double ref_moles, const double ref_temperature)
{
  ref_moles_           = ref_moles;
  inv_ref_moles_       = 1.0/ref_moles;
  ref_temperature_     = ref_temperature;
  inv_ref_temperature_ = 1.0/ref_temperature;
  if((ref_moles_ == 1.0) && (ref_temperature == 1.0)) {
    use_scaled_state_ = false;
  } else {
    use_scaled_state_ = true;
  }
}

void VariableVolumeReactor::Impl::GetReferenceScales(double *ref_moles, double *ref_temperature) const
{
  *ref_moles = ref_moles_;
  *ref_temperature = ref_temperature_;
}

double
  VariableVolumeReactor::Impl::GetMixtureSpecificHeat_Cv(const double state[])
{
  const int num_species            = GetNumSpecies();
  const double ref_temperature     = ref_temperature_;
  const double temperature         = ref_temperature*state[num_species+1];
  zerork::mechanism *mechanism_ptr;

  mechanism_ptr = GetMechanism();

  // the molar Cv calculation is normalized by the concentration sum. This 
  // means the mole fractions, or mole counts can also be used in place of
  // of the concentrations.  
  return mechanism_ptr->getMolarCvFromTC(temperature,
                                         &state[0],
                                         &specific_heats_[0]);
} 

// ---------------------------------------------------------------------------
// Public facing API
VariableVolumeReactor::VariableVolumeReactor(const char mechanism_name[],
                                             const char thermodynamics_name[],
                                             const char parser_log_name[],
                                             const MatrixType matrix_type,
                                             Volume *volume_function,
                                             HeatLoss *heat_loss_function)
{
  impl_ = new Impl(mechanism_name,
                   thermodynamics_name,
                   parser_log_name,
                   matrix_type,
                   volume_function,
                   heat_loss_function);             
}

VariableVolumeReactor::~VariableVolumeReactor()
{
  if(impl_ != NULL) {
    delete impl_;
  }
}
  
ReactorError
  VariableVolumeReactor::GetTimeDerivative(const double reactor_time,
                                           const double state[],
                                           double derivative[])
{
  return impl_->GetTimeDerivative(reactor_time,state,&derivative[0]);
}

ReactorError
  VariableVolumeReactor::GetChemicalHeatReleaseRate(const double reactor_time,
                                                    const double state[],
                                                    double *heat_release_rate)
{
  return impl_->GetChemicalHeatReleaseRate(reactor_time,
                                           state,
                                           heat_release_rate);
}

ReactorError VariableVolumeReactor::GetJacobian(const double reactor_time,
                                                const double state[],
                                                double jacobian[])
{
  return impl_->GetJacobian(reactor_time,state,&jacobian[0]);
}

int VariableVolumeReactor::GetJacobianPattern(int row_id[],
                                              int col_id[])
{
  return impl_->GetJacobianPattern(&row_id[0],&col_id[0]);
}

int VariableVolumeReactor::GetNumStates() const
{
  return impl_->GetNumStates();
}

int VariableVolumeReactor::GetNumSpecies() const
{
  return impl_->GetNumSpecies();
}

int VariableVolumeReactor::GetNumReactions() const
{
  return impl_->GetNumReactions();
}

int VariableVolumeReactor::GetNumSteps() const
{
  return impl_->GetNumSteps();
}

int VariableVolumeReactor::GetJacobianSize() const
{
  return impl_->GetJacobianSize();
}

MatrixType VariableVolumeReactor::GetMatrixType() const
{
  return impl_->GetMatrixType();
}

const char * VariableVolumeReactor::GetMechanismName() const
{
  return impl_->GetMechanismName();
}

const char * VariableVolumeReactor::GetThermodynamicsName() const
{
  return impl_->GetThermodynamicsName();
}

const char * VariableVolumeReactor::GetParserLogName() const
{
  return impl_->GetParserLogName();
}

const char * VariableVolumeReactor::GetReactorInfo() const
{
  return impl_->GetReactorInfo();
}

zerork::mechanism * VariableVolumeReactor::GetMechanism()
{
  return impl_->GetMechanism();
}

// species/reaction accessors
int VariableVolumeReactor::GetIdOfState(const char *state_name) const
{
  return impl_->GetIdOfState(state_name);
}

const char * VariableVolumeReactor::GetNameOfStateId(const int state_id) const
{
  return impl_->GetNameOfStateId(state_id);
}


// A-Factor sensitivity utilities
double VariableVolumeReactor::GetAMultiplierOfForwardReactionId(const int reaction_id) const
{
  return impl_->GetAMultiplierOfForwardReactionId(reaction_id);
}

ReactorError VariableVolumeReactor::SetAMultiplierOfForwardReactionId(const int reaction_id, const double a_multiplier)
{
  return impl_->SetAMultiplierOfForwardReactionId(reaction_id,a_multiplier);
}

double VariableVolumeReactor::GetAMultiplierOfReverseReactionId(const int reaction_id) const
{
  return impl_->GetAMultiplierOfReverseReactionId(reaction_id);
}

ReactorError VariableVolumeReactor::SetAMultiplierOfReverseReactionId(const int reaction_id, const double a_multiplier)
{
  return impl_->SetAMultiplierOfReverseReactionId(reaction_id, a_multiplier);
}

double VariableVolumeReactor::GetAMultiplierOfStepId(const int step_id) const
{
  return impl_->GetAMultiplierOfStepId(step_id);
}

ReactorError VariableVolumeReactor::SetAMultiplierOfStepId(const int step_id, const double a_multiplier)
{
  return impl_->SetAMultiplierOfStepId(step_id,a_multiplier);
}

void VariableVolumeReactor::SetReferenceScales(const double ref_moles,
                                               const double ref_temperature)
{
  impl_->SetReferenceScales(ref_moles,
                            ref_temperature);
}
void VariableVolumeReactor::GetReferenceScales(double *ref_moles,
                                               double *ref_temperature) const
{
  impl_->GetReferenceScales(ref_moles,
                            ref_temperature);
}

void VariableVolumeReactor::SetVolumeMultiplier(const double v_multiplier)
{
  impl_->SetVolumeMultiplier(v_multiplier);
}


double VariableVolumeReactor::GetVolumeMultiplier() const
{
  return impl_->GetVolumeMultiplier();
}

double
  VariableVolumeReactor::GetMixtureSpecificHeat_Cv(const double state[])
{
  return impl_->GetMixtureSpecificHeat_Cv(state);
}

void VariableVolumeReactor::GetSpeciesHydrogenCount(int num_atoms[]) const
{
  impl_->GetSpeciesHydrogenCount(num_atoms);
}
void VariableVolumeReactor::GetSpeciesNitrogenCount(int num_atoms[]) const
{
  impl_->GetSpeciesNitrogenCount(num_atoms);
}
void VariableVolumeReactor::GetSpeciesCarbonCount(int num_atoms[]) const
{
  impl_->GetSpeciesCarbonCount(num_atoms);
}
void VariableVolumeReactor::GetSpeciesOxygenCount(int num_atoms[]) const
{
  impl_->GetSpeciesOxygenCount(num_atoms);
}
void VariableVolumeReactor::GetSpeciesMolecularWeight(double molecular_weight[]) const
{
  impl_->GetSpeciesMolecularWeight(molecular_weight);
}

double VariableVolumeReactor::GetGasConstant() const
{
  return impl_->GetGasConstant();
}
