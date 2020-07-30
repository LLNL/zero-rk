#ifndef IDTSOLVERS_H
#define IDTSOLVERS_H

#include "idtControlParams.h"

int solveIdtSimple(idtControlParams *idtCtrl,
                   double *results,
                   double *solve_time);
int solveIdtOriginal(idtControlParams *idtCtrl,
                     double *results,
                     double *solveTime);

int solveIdtPerturbRxn(const int rxnId,
                       idtControlParams *idtCtrl,
                       double *results,
                       double *solveTime);


#endif
