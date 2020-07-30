#include <stdlib.h>
#include <stdio.h>

#include "event_counter.h"


typedef struct
{
  int event_id_;
  int count_;
} EventCount;

static int CompareCountDescending(const void *p1, const void *p2)
{
  EventCount *element1 = (EventCount *)p1;
  EventCount *element2 = (EventCount *)p2;
  if(element1->count_ > element2->count_) {
    return -1; // element1 goes before element2
  } else if (element1->count_ < element2->count_) {
    return  1; // element1 goes after element2
  }
  return 0;    // element1 and element2 are equivalent
}


EventCounter::EventCounter()
{
  count_map_.clear();
}

void EventCounter::AddEvent(const int event)
{
  std::map<int, int>::const_iterator iter;
  iter = count_map_.find(event);
  if(iter == count_map_.end()) {
    // new event
    count_map_[event] = 1;
  } else {
    ++count_map_[event];
  }
}

void EventCounter::GetEventList(std::vector<int> *events,  
                                std::vector<int> *counts) const
{
  events->clear();
  counts->clear();

  std::map<int, int>::const_iterator iter;
  for(iter = count_map_.begin(); iter != count_map_.end(); ++iter) {
    events->push_back(iter->first);
    counts->push_back(iter->second);
  }
}

void EventCounter::GetSortedEventList(std::vector<int> *events,  
                                      std::vector<int> *counts) const
{
  EventCount list_item;
  std::vector<EventCount> sort_list;
  GetEventList(events,counts);

  // Copy events and counts
  for(size_t j=0; j<events->size(); ++j) {
    list_item.event_id_ = events->at(j);
    list_item.count_    = counts->at(j);
    sort_list.push_back(list_item);
  }
  qsort(&sort_list[0],
        sort_list.size(),
        sizeof(EventCount),
        CompareCountDescending);
  // Copy sorted events and counts
  for(size_t j=0; j<events->size(); ++j) {
    events->at(j) = sort_list[j].event_id_;
    counts->at(j) = sort_list[j].count_;
  }      
    
}

int EventCounter::GetNumDistinctEvents() const
{
  return (int)count_map_.size();
}

int EventCounter::GetTotalCount() const
{
  int total_count = 0;
  std::map<int, int>::const_iterator iter;
  for(iter = count_map_.begin(); iter != count_map_.end(); ++iter) {
    total_count += iter->second;
  }
  return total_count;
}
 
