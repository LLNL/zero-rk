#include <math.h>

#include "atol_crossing.h"


AbsoluteToleranceCrossing::AbsoluteToleranceCrossing(
    const double absolute_tolerance)
{
  num_updates_ = 0;
  num_mass_fractions_ = 0;
  absolute_tolerance_ = absolute_tolerance;
  last_update_exceeds_.clear();
}

void AbsoluteToleranceCrossing::UpdateCurrentState(
    const CVReactorState &state)
{
  if(num_updates_ == 0) {
    // setup first state update
    num_mass_fractions_ = state.mass_fraction.size();
    last_update_exceeds_.assign(num_mass_fractions_, false);

    for(int j=0; j<num_mass_fractions_; ++j) {
      if(fabs(state.mass_fraction[j]) > absolute_tolerance_) {
        last_update_exceeds_[j] = true;
      }
    }
    ++num_updates_;
  } else {
    // check if any mass fraction has changed its past state
    if(num_mass_fractions_ == (int)state.mass_fraction.size()) {

      for(int j=0; j<num_mass_fractions_; ++j) {
        if(fabs(state.mass_fraction[j]) > absolute_tolerance_ && 
           last_update_exceeds_[j] == false) {

          last_update_exceeds_[j] = true;
          event_counter_.AddEvent(j);          
        } else if(fabs(state.mass_fraction[j]) <= absolute_tolerance_ && 
           last_update_exceeds_[j] == true) {

          last_update_exceeds_[j] = false;
          event_counter_.AddEvent(j);          
        }
      
      }
      ++num_updates_;
    } else {
      printf(
          "# WARNING: In AbsoluteToleranceCrossing::UpdateCurrentState(...),\n"
          "#          current state has %d species,\n"
          "#          but the class has %d species.\n"
          "#          Skipping update.\n",
          (int)state.mass_fraction.size(),
          num_mass_fractions_);
    }
  } // if(num_updates_ == 0) else
}


void AbsoluteToleranceCrossing::GetSortedCrossingList(
    std::vector<int> *species_id,
    std::vector<int> *counts) const
{
  species_id->clear();
  counts->clear();
  if(num_updates_ >= 2) {
    event_counter_.GetSortedEventList(species_id,counts);
  }
}
