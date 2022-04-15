
if(ENABLE_MPI)
set(COMMON_SRC    perturbAFactor_common.cpp matrix_funcs.cpp ode_funcs.cpp
                  cv_param_sparse.cpp utility_funcs.cpp GSA_AFactorIFP.cpp
                  idtSolvers.cpp idtControlParams.cpp gsa_stats.cpp)

add_mpi_executable(perturbAFactorGSA_mpi.x perturbAFactor_mpi.cpp ${COMMON_SRC})
target_link_libraries(perturbAFactorGSA_mpi.x zerorkutilities zerork
                      sundials_nvecserial sundials_cvode superlu spify)
install(TARGETS perturbAFactorGSA_mpi.x
        RUNTIME DESTINATION bin)
add_spifyIFP_target(GSA_AFactorIFP perturbAFactorGSA_mpi.x)
endif()

