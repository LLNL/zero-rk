#include <string>

#include <reactor/variable_volume_reactor.h>
#include <mechanism_info/mechanism_info.h>

#include "user_functions.h"

int VariableVolumeRHS(realtype t,        // [in] ODE system time
                      N_Vector y,        // [in] ODE state vector
                      N_Vector ydot,     // [out] ODE state derivative
		      void *params)      // [in/out]
{
  UserData *user_data = (UserData *)params;
  double *state      = NV_DATA_S(y);    // pointers to data array for N_Vector
  double *derivative = NV_DATA_S(ydot);
  const int num_states = user_data->GetReactor()->GetNumStates();
  // TODO: add error check
  user_data->GetReactor()->GetTimeDerivative(t,state,derivative);

 return 0;
}
int VariableVolumeRHS_Limit1(realtype t,        // [in] ODE system time
                      N_Vector y,        // [in] ODE state vector
                      N_Vector ydot,     // [out] ODE state derivative
		      void *params)      // [in/out]
{
  UserData *user_data = (UserData *)params;
  double *state      = NV_DATA_S(y);    // pointers to data array for N_Vector
  double *derivative = NV_DATA_S(ydot);
  const int num_species = user_data->GetReactor()->GetNumStates()-2;
  // TODO: add error check
  user_data->GetReactor()->GetTimeDerivative(t,state,derivative);
  for(int j=0; j<num_species; ++j) {
    //if(state[j] <= 1.0e-20 && derivative[j] < 0.0) {
    //  double x = state[j]*1.0e+20;
    //  double weight = (x < 0.0) ? 0.0 : x*x*(3.0-2.0*x);
    //  derivative[j] *= weight;
    //}
    if(state[j] <= 0.0 && derivative[j] < 0.0) {
      derivative[j] = 0.0;
    }
  }

  return 0;
}


int VariableVolumeRHS_Tmin(realtype t,        // [in] ODE system time
                      N_Vector y,        // [in] ODE state vector
                      N_Vector ydot,     // [out] ODE state derivative
		      void *params)      // [in/out]
{
  const double min_temperature = 600.0;
  UserData *user_data = (UserData *)params;
  double *state      = NV_DATA_S(y);    // pointers to data array for N_Vector
  double *derivative = NV_DATA_S(ydot);
  const int num_states = user_data->GetReactor()->GetNumStates();
  const double temperature = user_data->GetTemperature(state);
  // TODO: add error check
  user_data->GetReactor()->GetTimeDerivative(t,state,derivative);
  if(temperature < min_temperature) {
    // zero out any chemical change
    for(int j=0; j<num_states-1; ++j) {
      derivative[j] = 0.0;
    }
  }

  return 0;
}

// cvode version specific definition - TODO: remove the need to switch
// function definitions by hand
#if defined SUNDIALS2
int VariableVolumeDenseJacobian(long int N,     // [in] ODE system size
                                realtype t,     // [in] ODE system time
                                N_Vector y,     // [in] ODE state vector
                                N_Vector ydot,  // [in] ODE state derivative
                                DlsMat Jac,     // [out] ODE Jacobian
		                void *params,   // [in/out]
                                N_Vector tmp1,  // [out] N-length workspace
                                N_Vector tmp2,  // [out] N-length workspace
                                N_Vector tmp3)  // [out] N-length workspace
{
#elif defined SUNDIALS3 || defined SUNDIALS4
int VariableVolumeDenseJacobian(realtype t,     // [in] ODE system time
                                N_Vector y,     // [in] ODE state vector
                                N_Vector ydot,  // [in] ODE state derivative
                                SUNMatrix Jac,     // [out] ODE Jacobian
		                void *params,   // [in/out]
                                N_Vector tmp1,  // [out] N-length workspace
                                N_Vector tmp2,  // [out] N-length workspace
                                N_Vector tmp3)  // [out] N-length workspace
{
#endif
  UserData *user_data = (UserData *)params;
  double *state      = NV_DATA_S(y);    // pointers to data array for N_Vector
  //double *derivative = NV_DATA_S(ydot);

  const int num_nonzeros = user_data->GetReactor()->GetJacobianSize();
  std::vector<double> jacobian;
  std::vector<int> row_id, col_id;

  jacobian.assign(num_nonzeros,0.0);
  row_id.assign(num_nonzeros,0);
  col_id.assign(num_nonzeros,0);
  user_data->GetReactor()->GetJacobianPattern(&row_id[0],&col_id[0]);
  user_data->GetReactor()->GetJacobian(t,state,&jacobian[0]);

  // load non-zero elements, the Jacobian matrix is zeroed before entry
  for(int j=0; j<num_nonzeros; ++j) {
#ifdef SUNDIALS2
    DENSE_ELEM(Jac, row_id[j], col_id[j]) = jacobian[j];
#else
    SM_ELEMENT_D(Jac, row_id[j], col_id[j]) = jacobian[j];
#endif
  }
  return 0;
}

#if defined SUNDIALS2
int VariableVolumePreconditionerSetup(realtype t,// [in] ODE system time
                                      N_Vector y,      // [in] ODE state vector
                                      N_Vector ydot,   // [in] ODE state derivative
                                      booleantype jok,
                                      booleantype *new_j,
                                      realtype gamma,
                                      void *params,    // [in/out]
                                      N_Vector tmp1,   // [out] N-length workspace
                                      N_Vector tmp2,   // [out] N-length workspace
                                      N_Vector tmp3)   // [out] N-length workspace
{
#elif defined SUNDIALS3 || defined SUNDIALS4
  int VariableVolumePreconditionerSetup(realtype t,// [in] ODE system time
                                      N_Vector y,      // [in] ODE state vector
                                      N_Vector ydot,   // [in] ODE state derivative
                                      booleantype jok,
                                      booleantype *new_j,
                                      realtype gamma,
                                      void *params)    // [in/out]
{
#endif
  UserData *user_data    = (UserData *)params;
  const int num_nonzeros = user_data->GetReactor()->GetJacobianSize();
  const int num_states   = user_data->GetReactor()->GetNumStates();
  double *state          = NV_DATA_S(y); // pointers to data array for N_Vector
  int    *diagonal_id    = user_data->GetDiagonalId();
  double *jacobian       = user_data->GetJacobian();
  double *saved_jacobian = user_data->GetSavedJacobian();
  int error_flag;

  if(!jok) {
    // The Jacobian is not okay, need to recompute
    user_data->GetReactor()->GetJacobian(t,
                                         state,
                                         saved_jacobian);
    ++user_data->num_new_jacobians;
    (*new_j) = true;
  } else {
    (*new_j) = false;
  }

  // compute I - gamma*J
  // this could be updated with blas routines
  for(int j=0; j<num_nonzeros; ++j) {
    jacobian[j] = -gamma*saved_jacobian[j];
  }
  for(int j=0; j<num_states; ++j) {
    jacobian[diagonal_id[j]] += 1.0;
  }

  if(user_data->GetSparseMatrix()->IsFirstFactor()) {
    error_flag =
      user_data->GetSparseMatrix()->FactorNewPatternCCS(num_nonzeros,
						     user_data->GetRowId(),
						     user_data->GetColumnSum(),
                                                     jacobian);
  } else {
    error_flag =
      user_data->GetSparseMatrix()->FactorSamePattern(jacobian);
  }
  ++user_data->num_preconditioner_setups;
  return error_flag;
}

#if defined SUNDIALS2
int VariableVolumePreconditionerSolve(realtype t,// [in] ODE system time
                                N_Vector y,      // [in] ODE state vector
                                N_Vector ydot,   // [in] ODE state derivative
                                N_Vector r,      // [in] jacobian rhs
                                N_Vector z,      // [out]
				realtype gamma,
				realtype delta,
				int lr,
		                void *params,    // [in/out]
				N_Vector tmp)    // [out] N-length workspace
{
#elif defined SUNDIALS3 || defined SUNDIALS4
int VariableVolumePreconditionerSolve(realtype t,// [in] ODE system time
                                      N_Vector y,      // [in] ODE state vector
                                      N_Vector ydot,   // [in] ODE state derivative
                                      N_Vector r,      // [in] jacobian rhs
                                      N_Vector z,      // [out]
                                      realtype gamma,
                                      realtype delta,
                                      int lr,
                                      void *params)    // [in/out]
{
#endif
  UserData *user_data = (UserData *)params;
  double *rhs         = NV_DATA_S(r); // pointers to data array for N_Vector
  double *solution    = NV_DATA_S(z); // pointers to data array for N_Vector

  int error_flag = user_data->GetSparseMatrix()->Solve(rhs,solution);

  ++user_data->num_preconditioner_solves;

  return error_flag;
}
// ---------------------------------------------------------------------------
// UserData class member functions
UserData::UserData(const int task_num,
                   const char input_filename[],
                   const char volume_filename[],
                   const double initial_pressure,
                   const double initial_temperature,
                   const double fuel_fraction,
                   const CompositionMap &oxidizer_map)
{
  int error_flag;

  task_num_            = task_num;
  initial_pressure_    = initial_pressure;
  initial_temperature_ = initial_temperature;
  fuel_fraction_       = fuel_fraction;
  oxidizer_map_        = oxidizer_map;

  parser_ = new VariableVolumeBatchIFP(input_filename);
  if(parser_ == NULL) {
    printf("ERROR: In UserData constructor,\n");
    printf("       could not create input file parser from file %s\n",
           input_filename);
    fflush(stdout);
    exit(-1);
  }
  fuel_map_ = parser_->fuelComp();
  NormalizeCompositionMap(&fuel_map_);

  InterpolationType interpType = LINEAR_CLIPPED;
  if(parser_->volumeInterpolationMethod() == std::string("LINEAR_CLIPPED")) {
    interpType = LINEAR_CLIPPED;
  }
  if(parser_->volumeInterpolationMethod() == std::string("CUBIC_CLIPPED")) {
    interpType = CUBIC_CLIPPED;
  }

  heat_loss_ = new Adiabatic();
  volume_    = new VolumeFromFile(volume_filename,
                                  interpType,
                                  "!#", // comment characters to ignore
                                  parser_->volumeMult(),
                                  parser_->strokeMult(),
                                  true, // csv file
                                  1);   // skip first line

  reactor_ = new VariableVolumeReactor(parser_->mechFile().c_str(),
                                       parser_->thermFile().c_str(),
                                       parser_->mechLogFile().c_str(),
                                       //DENSE_COL_MAJOR,
                                       COMPRESSED_COL_STORAGE,
                                       volume_,
                                       heat_loss_);

  // TO DO: add checks to make sure all objects are created

  initial_time_ = volume_->GetMinTime();

  error_flag = BuildInitialState();
  if(error_flag != 0) {
    printf("ERROR: In UserData constructor,\n");
    printf("       could not BuildInitialState()\n");
    fflush(stdout);
    exit(-1);
  }

  int num_nonzeros = reactor_->GetJacobianPattern(NULL,NULL);

  reactor_->SetReferenceScales(reference_moles_,
                               parser_->refTemp());


  jacobian_.assign(num_nonzeros, 0.0);
  saved_jacobian_.assign(num_nonzeros, 0.0);

  row_id_.assign(num_nonzeros, 0);
  column_id_.assign(num_nonzeros, 0);

  column_sum_.assign(reactor_->GetNumStates()+1,0);
  diagonal_id_.assign(reactor_->GetNumStates(),-1); // set to negative one
                                                    // to check if we have
                                                    // a missing diagonal

  //printf("# Number of nonzeros in the Jacobian: %d\n",num_nonzeros);
  reactor_->GetJacobianPattern(&row_id_[0],&column_id_[0]);
  //for(int j=0; j<num_nonzeros; ++j) {
  //  printf("row %2d  col %2d\n",row_id_[j],column_id_[j]);
  //}
  //fflush(stdout);

  // assemble the column_sum_ array and the diagonal_id
  for(int j=0; j<num_nonzeros; ++j) {
    if(row_id_[j] == column_id_[j]) {
      diagonal_id_[row_id_[j]] = j;
    }
    ++column_sum_[column_id_[j]+1];
  }
  for(int j=0; j<reactor_->GetNumStates(); ++j) {
    column_sum_[j+1] += column_sum_[j];
  }
  // check to make sure we have the complete diagonal, since we need
  // to add an identity to the sparse matrix
  for(int j=0; j<reactor_->GetNumStates(); ++j) {
    if(diagonal_id_[j] == -1) {
      printf("ERROR: row %d does not have a diagonal term in the Jacobian\n",
             j);
      fflush(stdout);
      exit(-1);
    }
  }

  sparse_matrix_ = new SparseMatrix(reactor_->GetNumStates(),
                                    num_nonzeros);
  num_preconditioner_setups = 0;
  num_preconditioner_solves = 0;
  num_new_jacobians = 0;
}
UserData::~UserData()
{
  if(heat_loss_ != NULL) {
    delete heat_loss_;
  }
  if(volume_ != NULL) {
    delete volume_;
  }
  if(parser_ != NULL) {
    delete parser_;
  }
  if(reactor_ != NULL) {
    delete reactor_;
  }
  if(sparse_matrix_ != NULL) {
    delete sparse_matrix_;
  }
}

int UserData::BuildInitialState()
{
  const int num_species = reactor_->GetNumSpecies();
  const int num_states  = reactor_->GetNumStates();
  int species_id;
  std::string state_prefix("Moles_");
  std::string search_name;
  double initial_volume;
  double initial_moles;
  double dvolume_dt;
  double mole_frac_sum;

  std::map<std::string, double>::const_iterator map_iter;

  //printf("num_species: %d\n",num_species);
  //printf("num_states:  %d\n",num_states);
  //for(int j=0; j<num_states; ++j) {
  //  printf("State[%d]: %s\n",j,reactor_->GetNameOfStateId(j));
  //  fflush(stdout);
  //}

  //fflush(stdout);
  tracked_species_ids_.clear();
  initial_state_.assign(num_states,0.0);
  volume_->GetVolume(initial_time_,
                     NULL,
                     &initial_volume,
                     &dvolume_dt);

  initial_moles = (initial_pressure_*initial_volume)/
    (USER_GAS_CONSTANT*initial_temperature_);

  // parse the fuel composition information
  mole_frac_sum = 0.0;
  for(map_iter = fuel_map_.begin();
      map_iter != fuel_map_.end();
      ++map_iter) {

    search_name = map_iter->first;
    search_name.insert(0,state_prefix); // insert the prefix at the beginning
    species_id = reactor_->GetIdOfState(search_name.c_str());
    if(species_id < 0) {
      printf("ERROR: could not find state %s for fuelComp species %s\n",
             search_name.c_str(),map_iter->first.c_str());
      fflush(stdout);
      return -1;
    }
    tracked_species_ids_.push_back(species_id);
    initial_state_[species_id] = fuel_fraction_*map_iter->second;
    mole_frac_sum += fuel_fraction_*map_iter->second;
  }
  // parse the oxidizeer composition information
  for(map_iter = oxidizer_map_.begin();
      map_iter != oxidizer_map_.end();
      ++map_iter) {

    search_name = map_iter->first;
    search_name.insert(0,state_prefix); // insert the prefix at the beginning
    species_id = reactor_->GetIdOfState(search_name.c_str());
    if(species_id < 0) {
      printf("ERROR: could not find state %s for fuelComp species %s\n",
             search_name.c_str(),map_iter->first.c_str());
      fflush(stdout);
      return -1;
    }
    tracked_species_ids_.push_back(species_id);
    initial_state_[species_id] += (1.0-fuel_fraction_)*map_iter->second;
    mole_frac_sum += (1.0-fuel_fraction_)*map_iter->second;
  }

  // normalize the initial species compostion
  // may be MOLE_FRACTION or MASS_FRACTION
  for(int j=0; j<num_species; ++j) {
    initial_state_[j] *= 1.0/mole_frac_sum;
  }

  // the initial_state_ [0:num_species-1] is either mole or mass fractions
  // if mass fraction it need to be converted to mole fraction
  if(parser_->compType() == std::string("MASS_FRACTION")) {
    std::vector<double> converted_comp;
    MechanismInfo info(parser_->mechFile().c_str(),
                       parser_->thermFile().c_str(),
                       "");

    converted_comp.assign(num_species, 0.0);
    info.ConvertMassToMoleFraction(&initial_state_[0],
                                   &converted_comp[0]);

    printf("# Converting initial composition from mass to mole fractions:\n");
    for(int j=0; j<num_species; ++j) {
      initial_state_[j] = converted_comp[j];
    }
  } else if(parser_->compType() != std::string("MOLE_FRACTION")) {
    printf("WARNING: compType = %s not recognized\n",
           parser_->compType().c_str());
    printf("         assuming fuelComp is specified as mole fractions.\n");
    fflush(stdout);
  }


  // scale the total number of moles
  reference_moles_ = initial_moles; // parser_->refMoles();
  initial_moles *= (1.0/reference_moles_);
  for(int j=0; j<num_species; ++j) {
    initial_state_[j] *= initial_moles;
  }
  initial_state_[num_species]   = initial_moles;
  initial_state_[num_species+1] = initial_temperature_/parser_->refTemp();

  return 0;
}


void UserData::GetInitialState(double state[]) const
{
  const int num_states = reactor_->GetNumStates();
  for(int j=0; j<num_states; ++j) {
    state[j] = initial_state_[j];
  }
}

double UserData::GetVolume(const double t,
                           const double state[]) const
{
  double volume, dvolume_dt;
  volume_->GetVolume(t,
                     state,
                     &volume,
                     &dvolume_dt);
  return volume;
}

double UserData::GetVolumeRate(const double t,
                               const double state[]) const
{
  double volume, dvolume_dt;
  volume_->GetVolume(t,
                     state,
                     &volume,
                     &dvolume_dt);
  return dvolume_dt;
}

double UserData::GetTemperature(const double state[]) const
{
  const int num_species = reactor_->GetNumSpecies();
  double ref_moles, ref_temperature;
  reactor_->GetReferenceScales(&ref_moles,&ref_temperature);

  return state[num_species+1]*ref_temperature;
}


double UserData::GetTotalMoles(const double state[]) const
{
  const int num_species = reactor_->GetNumSpecies();
  double ref_moles, ref_temperature;
  reactor_->GetReferenceScales(&ref_moles,&ref_temperature);

  return state[num_species]*ref_moles;
}

double UserData::GetPressure(const double t,
                             const double state[]) const
{
  double volume = GetVolume(t,state);
  double temperature = GetTemperature(state);
  double total_moles = GetTotalMoles(state);

  return total_moles*USER_GAS_CONSTANT*temperature/volume;
}

void UserData::GetMoleFractions(const double state[],
                                double mole_fractions[]) const
{
  const int num_species = reactor_->GetNumSpecies();
  const double inv_total_moles = 1.0/state[num_species];
  // note that the species moles and total moles should have the same reference
  // scale, so they don't need to be rescaled for the mole fraction
  for(int j=0; j<num_species; ++j) {
    mole_fractions[j] = state[j]*inv_total_moles;
  }
}

double UserData::GetChemicalHeatReleaseRate(const double t,
                                            const double state[]) const
{
  double heat_release_rate;
  reactor_->GetChemicalHeatReleaseRate(t, state, &heat_release_rate);
  return heat_release_rate;
}

void UserData::GetTrackedSpeciesIds(std::vector<int> &species_ids) const
{
  const int num_tracked = tracked_species_ids_.size();
  species_ids.clear();
  for(int j=0; j<num_tracked; ++j) {
    species_ids.push_back(tracked_species_ids_[j]);
  }
}

void UserData::GetTrackedMoleFractions(const double state[],
                                     std::vector<double> &mole_fractions) const
{
  const int num_tracked = tracked_species_ids_.size();
  const int num_species = reactor_->GetNumSpecies();
  const double inv_total_moles = 1.0/state[num_species];
   mole_fractions.clear();
  for(int j=0; j<num_tracked; ++j) {
    mole_fractions.push_back(state[tracked_species_ids_[j]]*inv_total_moles);
  }
}

void NormalizeCompositionMap(CompositionMap *composition)
{
  double composition_sum = 0.0;
  CompositionMap::iterator map_iter;
  for(map_iter = composition->begin();
      map_iter != composition->end();
      ++map_iter) {
    //printf("map_iter->second = %.18g\n",map_iter->second); fflush(stdout);
    composition_sum += map_iter->second;
  }
  if(composition_sum != 0.0) {
    for(map_iter = composition->begin();
        map_iter != composition->end();
        ++map_iter) {
      map_iter->second /= composition_sum;
    }
  }
}
