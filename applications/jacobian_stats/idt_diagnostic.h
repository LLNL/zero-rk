#ifndef IDT_DIAGNOSTIC_H_
#define IDT_DIAGNOSTIC_H_

#include "event_counter.h"
#include "idtControlParams.h"

typedef struct
{
  double current_time_;
  double current_temperature_;
  double last_time_step_;
  int most_negative_species_id_;
  std::vector<int> error_species_id_;
} TimeHistoryElement;

void checkCommandLine(int inpArgc, char *inpArgv[]);

void getInitialActiveSpecies(BasicReactorIFP *parser,
                             zerork::mechanism *mech,
                             std::vector<int> *initial_active_species);

void writeMechanismReport(BasicReactorIFP *parser,
                          idtControlParams *idt_control);

int GetCVodeErrorSpecies(idtControlParams *idt_control,
                         const double threshold,
                         std::vector<int> *species_id,
                         std::vector<double> *species_error);

int GetCVodeErrorReactions(idtControlParams *idt_control,
                           const double threshold,
                           std::vector<int> *reaction_id,
                           std::vector<double> *reaction_error);

int GetErrorHistoryRanges(const std::vector<TimeHistoryElement> &error_history,
                          const int species_id,
                          double *min_time, 
                          double *max_time,
                          double *min_temperature, 
                          double *max_temperature,
                          double *min_time_step,
                          double *max_time_step,
                          double *mean_time_step);  
       



#endif
