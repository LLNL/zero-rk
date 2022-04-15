
set(COMMON_SRCS VariableVolumeBatchIFP.cpp
                complete_solver.cpp
                volume_functions.cpp sparse_matrix.cpp
                user_functions.cpp)


add_executable(variable_volume_batch.x variable_volume_batch.cpp ${COMMON_SRCS})
target_link_libraries(variable_volume_batch.x zerorkutilities mechanisminfo reactor zerork
                      sundials_nvecserial sundials_cvode superlu spify)
install(TARGETS variable_volume_batch.x
        RUNTIME DESTINATION bin)
set(DRIVERS variable_volume_batch.x)

if(ENABLE_MPI)
add_mpi_executable(variable_volume_batch_mpi.x variable_volume_batch_mpi.cpp ${COMMON_SRCS})
target_link_libraries(variable_volume_batch_mpi.x zerorkutilities mechanisminfo reactor zerork
                      sundials_nvecserial sundials_cvode superlu spify)
install(TARGETS variable_volume_batch_mpi.x
        RUNTIME DESTINATION bin)
set(DRIVERS ${DRIVERS} variable_volume_batch_mpi.x)
endif()

add_spifyIFP_target(VariableVolumeBatchIFP ${DRIVERS})

