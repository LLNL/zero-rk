
set(MAIN_SRC      cvIDT_sweep.cpp)
set(COMMON_SRC    matrix_funcs.cpp ode_funcs.cpp cv_param_sparse.cpp utility_funcs.cpp)
set(SPIFY_SRC     sweep_util_yml.cpp idt_sweep_IFP.cpp)
set(SPIFY_APPS    constVolumeWSR.x)

add_executable(constVolumeWSR.x ${MAIN_SRC} ${COMMON_SRC} ${SPIFY_SRC})
target_link_libraries(constVolumeWSR.x zerork sundials_nvecserial sundials_cvode superlu spify)
target_compile_definitions(constVolumeWSR.x PRIVATE SPIFY)
install(TARGETS constVolumeWSR.x
        RUNTIME DESTINATION bin)

if(ENABLE_MPI)
add_mpi_executable(constVolumeWSR_mpi.x ${MAIN_SRC} ${COMMON_SRC} ${SPIFY_SRC})
target_link_libraries(constVolumeWSR_mpi.x zerork sundials_nvecserial sundials_cvode superlu spify)
target_compile_definitions(constVolumeWSR_mpi.x PRIVATE SPIFY)
install(TARGETS constVolumeWSR_mpi.x
        RUNTIME DESTINATION bin)
set(SPIFY_APPS    "${SPIFY_APPS};constVolumeWSR_mpi.x")
endif()

add_spifyIFP_target(idt_sweep_IFP ${SPIFY_APPS})
