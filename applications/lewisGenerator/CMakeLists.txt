
if(ENABLE_MPI)
add_mpi_executable(lewisGenerator.x lewisGenerator.cpp
                   compute_lewis.cpp set_initial_conditions.cpp
                   flame_params.cpp UnsteadyFlameIFP.cpp)

target_link_libraries(lewisGenerator.x zerorkutilities mechanisminfo
                      reactor zerorktransport zerork spify sundials_nvecserial)
target_link_libraries(lewisGenerator.x zerork_equilibrium)

install(TARGETS lewisGenerator.x
        RUNTIME DESTINATION bin)
add_spifyIFP_target(UnsteadyFlameIFP lewisGenerator.x)
endif()
