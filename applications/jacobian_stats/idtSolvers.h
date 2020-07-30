#ifndef IDTSOLVERS_H
#define IDTSOLVERS_H

#include "idtControlParams.h"

int solveIdtSimple(idtControlParams *idtCtrl,
                   double idt[],
                   double *solveTime);

int solveIdtOriginal(idtControlParams *idtCtrl,
                     double idt[],
                     double *solveTime);

int solveIdtPerturbRxn(const int rxnId,
                       idtControlParams *idtCtrl,
                       double idt[],
                       double *solveTime);

typedef struct
{
  double pressure;
  double temperature;
  double density;
  double molecular_weight;

  std::vector<double> mole_fraction;
  std::vector<double> mass_fraction;
  std::vector<double> concentration;
} CVReactorState;


int marchIdtStep(idtControlParams *idtCtrl,
                 double t_start,
                 double *t_stop,
                 double *solveTime);

int marchIdtOneStep(idtControlParams *idtCtrl,
                    double t_start,
                    double *t_stop,
                    double *solveTime);

void checkAllFlagsFromIdtStep(idtControlParams *idtCtrl,
                              double t_current,
                              int cvode_flag,
                              double idt[],
                              FILE *check_file);

int checkNextIdtStep(idtControlParams *idtCtrl,
                     double t_current,
                     int verbosity);

void initializeIdtStep(idtControlParams *idtCtrl,
                      double idt[]);

void getState(idtControlParams *idtCtrl,
              CVReactorState *state);

void getDerivativesFromState(idtControlParams *idtCtrl,
                             CVReactorState *state,
                             std::vector<double> &net_species_rate,
                             std::vector<double> &step_rate);

void writeCVodeWarningsErrors(idtControlParams *idtCtrl,
                              double t_current,
                              int cvode_flag,
                              FILE *check_file);

void writeCVodeReport(idtControlParams *idtCtrl,
                    int cvode_flag,
                    FILE *fptr);
void writeThermoStateReport(idtControlParams *idtCtrl,
                            double t_current,
                            FILE *fptr);
void writeOdeStateReport(idtControlParams *idtCtrl,
                         double t_current,
                         FILE *fptr);

void printfDashedLine(); 
void fprintfDashedLine(FILE *fptr);

double GetMinMoleFraction(const CVReactorState &state,
                          int *min_species_id);
double GetNegativeMoleFractionSum(const CVReactorState &state);




#endif
