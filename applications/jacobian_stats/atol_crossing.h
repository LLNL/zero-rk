#ifndef ATOL_CROSSING_H_
#define ATOL_CROSSING_H_

#include "idtSolvers.h"    // class CVReactorState
#include "event_counter.h" // class EventCounter

class AbsoluteToleranceCrossing
{
 public:
  explicit AbsoluteToleranceCrossing(const double absolute_tolerance);
  ~AbsoluteToleranceCrossing() {};

  void UpdateCurrentState(const CVReactorState &state);
  void GetSortedCrossingList(std::vector<int> *species_id,
                             std::vector<int> *counts) const;
  int GetNumUpdates() const {return num_updates_;}

 private:
  int num_updates_;
  int num_mass_fractions_;
  double absolute_tolerance_;
  std::vector<bool> last_update_exceeds_;
  EventCounter event_counter_;
};



#endif
