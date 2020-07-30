
set(MAIN_SRC      cvIDT_sweep.cpp)
set(MAIN_MPI_SRC  cvIDT_sweep_mpi.cpp)
set(COMMON_SRC    matrix_funcs.cpp ode_funcs.cpp cv_param_sparse.cpp utility_funcs.cpp)
set(FORMATTED_SRC sweep_util_formatted.cpp)
set(SPIFY_SRC     sweep_util_yml.cpp idt_sweep_IFP.cpp)
set(SPIFY_APPS    constVolumeWSR_yml.x)

add_executable(constVolumeWSR.x ${MAIN_SRC} ${COMMON_SRC} ${FORMATTED_SRC})
target_link_libraries(constVolumeWSR.x zerork sundials_nvecserial sundials_cvode superlu)
install(TARGETS constVolumeWSR.x)

add_executable(constVolumeWSR_yml.x ${MAIN_SRC} ${COMMON_SRC} ${SPIFY_SRC})
target_link_libraries(constVolumeWSR_yml.x zerork sundials_nvecserial sundials_cvode superlu spify)
target_compile_definitions(constVolumeWSR_yml.x PRIVATE SPIFY)
install(TARGETS constVolumeWSR_yml.x)

if(ENABLE_MPI)
add_mpi_executable(constVolumeWSR_mpi.x ${MAIN_MPI_SRC} ${COMMON_SRC} ${FORMATTED_SRC})
target_link_libraries(constVolumeWSR_mpi.x zerork sundials_nvecserial sundials_cvode superlu)
install(TARGETS constVolumeWSR_mpi.x)

add_mpi_executable(constVolumeWSR_mpi_yml.x ${MAIN_MPI_SRC} ${COMMON_SRC} ${SPIFY_SRC})
target_link_libraries(constVolumeWSR_mpi_yml.x zerork sundials_nvecserial sundials_cvode superlu spify)
target_compile_definitions(constVolumeWSR_mpi_yml.x PRIVATE SPIFY)
install(TARGETS constVolumeWSR_mpi_yml.x)
set(SPIFY_APPS    "${SPIFY_APPS};constVolumeWSR_mpi_yml.x")
endif()

add_spifyIFP_target(idt_sweep_IFP ${SPIFY_APPS})
