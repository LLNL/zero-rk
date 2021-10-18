#ifndef IDTSOLVERS_H
#define IDTSOLVERS_H

#include "idtControlParams.h"

int solveIdtSimple(idtControlParams *idtCtrl,
                   double idt[],
                   double *solveTime);
int solveIdtOriginal(idtControlParams *idtCtrl,
                     double idt[],
                     double *solveTime);

// int solveIdtPerturbRxn(const int rxnId,
//                        idtControlParams *idtCtrl,
//                        double idt[],
//                        double *solveTime);
int solveIdtGsaPerturbRxn(const double afactor_mult[],
                          idtControlParams *idtCtrl,
                          double idt[],
                          double *solveTime);


#endif
