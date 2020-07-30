#ifndef PERTURBAFACTOR_COMMON_H
#define PERTURBAFACTOR_COMMON_H

#include <string>

#include "idtControlParams.h"

void checkCommandLine(int inpArgc, char *inpArgv[]);
void getHeaderInfo(int inpArgc,
                   char *inpArgv[],
                   idtControlParams *ctrl,
                   bool doFullReport,
                   string &header);
void getColHeaderInfo(idtControlParams *ctrl,
                      string &header);
void getColHeaderInfo_sensitivity(idtControlParams *ctrl,
                                  string &header);

void calcRxnSensitivity(const idtControlParams *ctrl,
                        const double idtOrig[],
                        const double idtPerturb[],
                        double rxnSens[],
                        int sortedRxnIdx[]);

typedef struct
{
  int rxnId;
  double maxRelSens;
} maxRxnSens_t;
int compare_maxRxnSens_t(const void *A, const void *B);

#endif
