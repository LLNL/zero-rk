
if(ENABLE_MPI)
set(SRCS VariableVolumeGSAIFP.cpp
         variable_volume_gsa.cpp
         complete_solver_gsa.cpp
         volume_functions.cpp sparse_matrix.cpp
         user_functions_gsa.cpp
         gsa_stats.cpp integrated_function.cpp)

add_mpi_executable(variable_volume_gsa.x ${SRCS})
target_link_libraries(variable_volume_gsa.x mechanisminfo reactor zerork
                      zerorkutilities
                      sundials_nvecserial sundials_cvodes
                      superlu spify)
install(TARGETS variable_volume_gsa.x
        RUNTIME DESTINATION bin)

add_spifyIFP_target(VariableVolumeGSAIFP variable_volume_gsa.x)
endif()
