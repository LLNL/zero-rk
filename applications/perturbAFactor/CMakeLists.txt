
set(COMMON_SRC    perturbAFactor_common.cpp matrix_funcs.cpp ode_funcs.cpp
                  cv_param_sparse.cpp utility_funcs.cpp AFactorIFP.cpp
                  idtSolvers.cpp idtControlParams.cpp)
set(SPIFY_APPS    perturbAFactor_serial.x)

add_executable(perturbAFactor_serial.x perturbAFactor_serial.cpp ${COMMON_SRC})
target_link_libraries(perturbAFactor_serial.x zerork sundials_nvecserial sundials_cvode superlu spify)
install(TARGETS perturbAFactor_serial.x
        RUNTIME DESTINATION bin)

if(ENABLE_MPI)
add_mpi_executable(perturbAFactor_mpi.x perturbAFactor_mpi.cpp ${COMMON_SRC})
target_link_libraries(perturbAFactor_mpi.x zerork sundials_nvecserial sundials_cvode superlu spify)
install(TARGETS perturbAFactor_mpi.x
        RUNTIME DESTINATION bin)
set(SPIFY_APPS    "${SPIFY_APPS};perturbAFactor_mpi.x")
endif()

add_spifyIFP_target(AFactorIFP ${SPIFY_APPS})
