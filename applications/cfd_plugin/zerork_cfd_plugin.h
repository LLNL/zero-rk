#ifndef ZERORK_CVREACTOR_LIB_H
#define ZERORK_CVREACTOR_LIB_H


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

void zerork_cfd_plugin_setup_defaults(double rtol, double atol,
                       double sparse_threshold);

void zerork_cfd_plugin_close(void);


void zerork_cfd_plugin_setup_full
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
    const char* mech_filename,
    const char* therm_filename,
    const char* ck_log_filename,
    const char* reactor_log_filename,
    int* multireac,
    int* gpu_id
);

void zerork_solve_reactors(const int nReactors, const double dt,
                        const double *T, const double *P,
                        double *massFracPtr, double* cost, double* gpu);

long int zerork_cfd_plugin_solve(const double dt, const double T, const double P, double *massFracPtr);

void zerork_print_stats();
void zerork_reset_stats();

#ifdef __cplusplus
}
#endif


#endif
