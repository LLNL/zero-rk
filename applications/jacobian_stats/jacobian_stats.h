#ifndef JACOBIAN_STATS_H_
#define JACOBIAN_STATS_H_

#include "utilities/distribution.h"
#include "idtControlParams.h"

const int temperature_id = -1;

typedef struct
{
  int species_id;
  int step_id;
  int elementary_step_order;
  bool is_pressure_dependent;
  double value;
} ROPDerivative;

int compareROPDerivative(const void *aptr, const void *bptr);
int compareROPDerivativeDescending(const void *aptr, const void *bptr);

int initializeROPDerivativeVector(idtControlParams *idtCtrl,
                                  std::vector<ROPDerivative> *vec);
double updateROPDerivativeVector(const double rate_of_progress[],
                                 const double concentration[],
                                 std::vector<ROPDerivative> *vec);
void sortROPDerivativeVector(std::vector<ROPDerivative> *vec);

void getROPDerivativeVector(idtControlParams *idtCtrl,
                            std::vector<ROPDerivative> *vec,
                            double *min_concentration);

int getROPDerivative(idtControlParams *idtCtrl,
                     double *min_concentration,
                     double rate_of_progress[],
                     int species_id[],
                     int step_id[]);

void writeROPDerivativeReport(idtControlParams *idtCtrl,
                              const double t_current,
                              const double min_value,
                              FILE *fptr);

void getROPTempDerivativeVector(idtControlParams *idtCtrl,
                               std::vector<ROPDerivative> *vec,
                               double *min_concentration);

int getROPTempDerivative(idtControlParams *idtCtrl,
                         double *min_concentration,
                         double rate_of_progress[],
                         int step_id[]);

void writeROPTempDerivativeReport(idtControlParams *idtCtrl,
                                  const double t_current,
                                  const double min_value,
                                  FILE *fptr);

void updateROPDerivativeDistribution(idtControlParams *idtCtrl,
                                     const double t_current,
                                     zerork::utilities::Distribution *order1_cum,
                                     zerork::utilities::Distribution *order2_cum,
                                     zerork::utilities::Distribution *order3g_cum,
                                     zerork::utilities::Distribution *thirdfall_cum,
                                     zerork::utilities::Distribution *temp_cum,
                                     FILE *fptr);


#endif
