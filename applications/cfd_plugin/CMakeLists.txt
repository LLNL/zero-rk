

set(COMMON_SRCS matrix_funcs.cpp ode_funcs.cpp 
       cv_param_sparse.cpp utility_funcs.cpp)

if(ENABLE_MPI)
add_mpi_library(zerork_cfd_plugin ${COMMON_SRCS} zerork_cfd_plugin.cpp)
target_compile_definitions(zerork_cfd_plugin PRIVATE USE_MPI)
add_mpi_executable(zerork_cfd_plugin_tester.x zerork_cfd_plugin_tester.cpp ZeroRKCFDPluginTesterIFP.cpp)
target_compile_definitions(zerork_cfd_plugin_tester.x PRIVATE USE_MPI)
else()
add_library(zerork_cfd_plugin ${COMMON_SRCS} zerork_cfd_plugin.cpp)
add_executable(zerork_cfd_plugin_tester.x zerork_cfd_plugin_tester.cpp ZeroRKCFDPluginTesterIFP.cpp)
endif()

add_spifyIFP_target(ZeroRKCFDPluginTesterIFP zerork_cfd_plugin_tester.x)

target_link_libraries(zerork_cfd_plugin zerork superlu spify sundials_cvode sundials_nvecserial)
target_link_libraries(zerork_cfd_plugin_tester.x zerork_cfd_plugin)

install(TARGETS zerork_cfd_plugin zerork_cfd_plugin_tester.x)
