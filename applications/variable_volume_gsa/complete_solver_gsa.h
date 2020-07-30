#ifndef COMPLETE_SOLVER_H_
#define COMPLETE_SOLVER_H_

#include "user_functions_gsa.h"


const int NUM_FIXED_RESULTS = 6;
//  results[0] = current_time;
//  results[1] = max_dp_dt;
//  results[2] = time_max_dp_dt;
//  results[3] = total_heat_release;
//  results[4] = max heat release rate (piece-wise linear)
//  results[5] = max heat release rate (piece-wise backward quadratic)

typedef struct
{
  bool record_history;
  bool echo_stdout;
  int step_print_period;
  double min_time_step;

} TimeHistoryParams;

typedef struct aux_cvode_structs_t
{
  aux_cvode_structs_t() : alloc(0) {}
  int alloc;
#if defined SUNDIALS3 || defined SUNDIALS4
  SUNMatrix A;
  SUNLinearSolver LS;
#endif
#if defined SUNDIALS4
  SUNNonlinearSolver NLS;
#endif
} aux_cvode_structs_t;


void SetupCompleteSolver(void *cvode_memory,
                         aux_cvode_structs_t &aux_cvode,
                         UserData *user_data,
                         N_Vector cvode_state);
void SetupSparseSolver(void *cvode_memory,
                       aux_cvode_structs_t &aux_cvode,
                       UserData *user_data,
                       N_Vector cvode_state);
void SetupDenseSolver(void *cvode_memory,
                      aux_cvode_structs_t &aux_cvode,
                      UserData *user_data,
                      N_Vector cvode_state);

void DestroyAuxCVodeStructs(aux_cvode_structs_t &aux_cvode);

int SolveVariableVolume(void *cvode_memory,
                        UserData *user_data,
                        N_Vector cvode_state,
                        const double a_multipliers[],
                        const TimeHistoryParams &thist_params,
                        std::string *thist,
                        double *results,
                        double *simulation_time);

void BuildTimeHistoryLine(void *cvode_memory,
                          const double t,
                          N_Vector state,
                          UserData *user_data,
                          std::string &line);

void BuildTimeHistoryHeader(void *cvode_memory,
                            const double t,
                            N_Vector state,
                            UserData *user_data,
                            std::string &header);

void PrintStateDerivatives(const double t,
                           N_Vector state,
                           UserData *user_data);

int CheckFlag(void *flagvalue, const char *funcname, int opt);

double GetHighResolutionTime(void);
#endif
