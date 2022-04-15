
if(ENABLE_MPI)
add_mpi_executable(premixed_steady_flame_solver.x premixed_steady_flame_solver.cpp
                   kinsol_functions.cpp set_initial_conditions.cpp flame_params.cpp sparse_matrix.cpp
                   sparse_matrix_dist.cpp UnsteadyFlameIFP.cpp)

target_link_libraries(premixed_steady_flame_solver.x zerorkutilities mechanisminfo reactor
                      zerorktransport zerork superlu_dist superlu spify sundials_kinsol sundials_nvecparallel)

install(TARGETS premixed_steady_flame_solver.x
        RUNTIME DESTINATION bin)

add_spifyIFP_target(UnsteadyFlameIFP premixed_steady_flame_solver.x)
endif()
