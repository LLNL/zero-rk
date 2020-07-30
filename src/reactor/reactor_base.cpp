#include "reactor_base.h"


ReactorError ReactorBase::BuildMechanism(const char mechanism_name[],
                                         const char thermodynamics_name[],
                                         const char parser_log_name[])
{
  mechanism_name_      = std::string(mechanism_name);
  thermodynamics_name_ = std::string(thermodynamics_name);
  parser_log_name_     = std::string(parser_log_name);

  // TODO: add more robust file/validity checks
  mechanism_ = new zerork::mechanism(mechanism_name,
                                     thermodynamics_name,
                                     parser_log_name);
  if(mechanism_ == NULL) {
    return INVALID_MECHANISM;
  }
  num_species_ = mechanism_->getNumSpecies();
  if(num_species_ < 1) {
    return INVALID_MECHANISM;
  }
  num_reactions_ = mechanism_->getNumReactions();
  if(num_reactions_ < 1) {
    return INVALID_MECHANISM;
  }
  num_steps_ = mechanism_->getNumSteps();
  if(num_steps_ < 1) {
    return INVALID_MECHANISM;
  }
  // Initialize all A-Factors to one
  a_multipliers_.assign(num_steps_,1.0);
  return NONE;
}

void ReactorBase::DestroyMechanism()
{
  if(mechanism_ != NULL) {
    delete mechanism_;
  }
}

ReactorError ReactorBase::SetJacobianSize(const int jacobian_size)
{
  if(jacobian_size < 1) {
    return INDEX_OUT_OF_RANGE;
  }
  jacobian_size_ = jacobian_size;
  return NONE;
}



// Needs to set the following private data members
//  num_states_
//  state_names_
//  state_names_map_
int
  ReactorBase::BuildStateNamesMap(const std::vector<std::string> &state_names)
{
  num_states_ = static_cast<int>(state_names.size());
  //printf("num_states_: %d\n",num_states_);
  fflush(stdout);
  state_names_.clear();
  state_names_map_.clear();
  for(int j=0; j<num_states_; ++j) {
    state_names_.push_back(state_names[j]);
    state_names_map_[state_names[j]] = j;
  }
  return num_states_;
}


int ReactorBase::GetIdOfState(const char *state_name) const
{
  std::string search_name = std::string(state_name);
  std::map<std::string,int>::const_iterator iter;

  iter = state_names_map_.find(search_name);
  if(iter == state_names_map_.end()) {
    // state_name not found
    return -1;
  }
  return iter->second;
}
const char * ReactorBase::GetNameOfStateId(const int state_id) const
{
  if(state_id < 0 || state_id >= num_states_) {
    return "invalid state index";
  }
  return state_names_[state_id].c_str();
}



double 
  ReactorBase::GetAMultiplierOfForwardReactionId(const int reaction_id) const
{
  if(reaction_id < 0 || reaction_id >= num_reactions_) {
    // reaction index does not exist
    return 0.0;
  }
  const int step_id = mechanism_->getStepIdxOfRxn(reaction_id,1);
  if(step_id < 0 || step_id >= num_steps_) {
    // step index does not exist
    return 0.0;
  }
  return a_multipliers_[step_id];
}
ReactorError 
  ReactorBase::SetAMultiplierOfForwardReactionId(const int reaction_id, 
                                                 const double a_multiplier)
{
  if(reaction_id < 0 || reaction_id >= num_reactions_) {
    // reaction index does not exist to be set
    return INDEX_OUT_OF_RANGE;
  }
  const int step_id = mechanism_->getStepIdxOfRxn(reaction_id,1);
  if(step_id < 0 || step_id >= num_steps_) {
    // step index does not exist to be set
    return INDEX_OUT_OF_RANGE;
  }
  a_multipliers_[step_id] = a_multiplier;
  return NONE;
}

double 
  ReactorBase::GetAMultiplierOfReverseReactionId(const int reaction_id) const
{
  if(reaction_id < 0 || reaction_id >= num_reactions_) {
    // reaction index does not exist
    return 0.0;
  }
  const int step_id = mechanism_->getStepIdxOfRxn(reaction_id,-1);
  if(step_id < 0 || step_id >= num_steps_) {
    // step index does not exist
    return 0.0;
  }
  return a_multipliers_[step_id];
}

ReactorError 
  ReactorBase::SetAMultiplierOfReverseReactionId(const int reaction_id, 
                                                 const double a_multiplier)
{
  if(reaction_id < 0 || reaction_id >= num_reactions_) {
    // reaction index does not exist to be set
    return INDEX_OUT_OF_RANGE;
  }
  const int step_id = mechanism_->getStepIdxOfRxn(reaction_id,-1);
  if(step_id < 0 || step_id >= num_steps_) {
    // step index does not exist to be set
    return INDEX_OUT_OF_RANGE;
  }
  a_multipliers_[step_id] = a_multiplier;
  return NONE;
}

double 
  ReactorBase::GetAMultiplierOfStepId(const int step_id) const
{
  if(step_id < 0 || step_id >= num_steps_) {
    // step index does not exist
    return 0.0;
  }
  return a_multipliers_[step_id];
}
ReactorError 
  ReactorBase::SetAMultiplierOfStepId(const int step_id, 
                                      const double a_multiplier)
{
  if(step_id < 0 || step_id >= num_steps_) {
    // step index does not exist to be set
    return INDEX_OUT_OF_RANGE;
  }
  a_multipliers_[step_id] = a_multiplier;
  return NONE;
}
