#ifndef EVENT_COUNTER_H_
#define EVENT_COUNTER_H_

#include <vector>
#include <map>

// TODO: template by event type
class EventCounter{

 public:
  EventCounter();
  // no allocated data

  void AddEvent(const int event);
  void GetEventList(std::vector<int> *events,  
                    std::vector<int> *counts) const;
  // the event list is sorted by highest counts to lowest
  void GetSortedEventList(std::vector<int> *events,  
                          std::vector<int> *counts) const;
  int GetNumDistinctEvents() const;
  int GetTotalCount() const;
 private:
  std::map<int, int> count_map_;

};

#endif

