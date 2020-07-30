
if(ENABLE_MPI)
add_mpi_executable(premixed_unsteady_flame_solver.x premixed_unsteady_flame_solver.cpp
                   cvode_functions.cpp set_initial_conditions.cpp
                   flame_params.cpp sparse_matrix.cpp UnsteadyFlameIFP.cpp)

target_link_libraries(premixed_unsteady_flame_solver.x zerorkutilities mechanisminfo reactor
                      zerorktransport zerork superlu spify sundials_cvode sundials_nvecparallel)

install(TARGETS premixed_unsteady_flame_solver.x)

add_spifyIFP_target(UnsteadyFlameIFP premixed_unsteady_flame_solver.x)
endif()
