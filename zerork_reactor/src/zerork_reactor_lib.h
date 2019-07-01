#ifndef ZERORK_CVREACTOR_LIB_H
#define ZERORK_CVREACTOR_LIB_H


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


void zerork_reactor_setup_defaults(double rtol, double atol,
                       double sparse_threshold);

void zerork_reactor_close(void);


void zerork_reactor_setup_full
(
    int verbosity,
    int maxsteps,
    int maxord,
    double rtol,
    double atol,
    int abstol_dens,
    double maxdt_internal,
    double sparsethresh,
    double nlconvcoef,
    double epslin,
    int lsolver,
    int CVCP,
    const char* mech_filename,
    const char* therm_filename,
    const char* ck_log_filename,
    const char* reactor_log_filename,
    int n_reactors_max,
    int n_reactors_min,
    double logiA, double logiK, double logiQ,
    double logiM, double logiB, double logirNu,
    int glu_factor_alg, int glu_solve_alg,
    int* multireac,
    int* gpu_id
);

void zerork_solve_reactors_mz(const int nReactors, const double dt,
                           const double *T, const double *P,
                           const double* dpdt, const double* volume,
                           double *massFracPtr, double* cost, double* gpu);

void zerork_solve_reactors(const int nReactors, const double dt,
                        const double *T, const double *P,
                        const double* dpdt,
                        double *massFracPtr, double* cost, double* gpu);

long int zerork_reactor_solve(const double dt, const double T, const double P,
                           const double dpdt, double *massFracPtr);

long int zerork_reactor_solve_mr(const int nReactors, const double dt,
                              const double *T, const double *P,
                              const double* dpdt, double *massFracPtr,
                              double* reactorCost);


void zerork_print_stats();
void zerork_reset_stats();


#ifdef __cplusplus
}
#endif


#endif
