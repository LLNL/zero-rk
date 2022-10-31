


add_executable(idt_diagnostic.x BasicReactorIFP.cpp
               sparse_eigenvalues.cpp event_counter.cpp atol_crossing.cpp
               utility_funcs.cpp ode_funcs.cpp matrix_funcs.cpp cv_param_sparse.cpp
               idtControlParams.cpp idtSolvers.cpp special_species.cpp
               mechanism_stats.cpp jacobian_stats.cpp idt_diagnostic.cpp)

target_link_libraries(idt_diagnostic.x zerork zerorkutilities zerorktransport
                      sundials_nvecserial sundials_cvode superlu spify)
add_spifyIFP_target(BasicReactorIFP idt_diagnostic.x)

install(TARGETS idt_diagnostic.x
        RUNTIME DESTINATION bin)

