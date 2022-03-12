

set(COMMON_SRCS utility_funcs.cpp
       optionable.cpp reactor_base.cpp
       reactor_constant_volume.cpp
       reactor_constant_pressure.cpp
       reactor_nvector_serial.cpp solver_cvode.cpp 
       solver_seulex.cpp utility_funcs.cpp
       zerork_reactor_manager.cpp
       interfaces/superlu_manager/superlu_manager.cpp
       interfaces/lapack_manager/lapack_manager.cpp 
       interfaces/seulex_cpp/seulex_cpp.cpp)

if(ENABLE_MPI)
add_mpi_library(zerork_cfd_plugin SHARED ${COMMON_SRCS} zerork_cfd_plugin.cpp)
target_compile_definitions(zerork_cfd_plugin PRIVATE USE_MPI)
add_mpi_executable(zerork_cfd_plugin_tester.x zerork_cfd_plugin_tester.cpp ZeroRKCFDPluginTesterIFP.cpp)
target_compile_definitions(zerork_cfd_plugin_tester.x PRIVATE USE_MPI)
else()
add_library(zerork_cfd_plugin SHARED ${COMMON_SRCS} zerork_cfd_plugin.cpp)
add_executable(zerork_cfd_plugin_tester.x zerork_cfd_plugin_tester.cpp ZeroRKCFDPluginTesterIFP.cpp)
endif()

target_link_libraries(zerork_cfd_plugin zerork superlu spify sundials_cvode sundials_nvecserial)
target_link_options(zerork_cfd_plugin PRIVATE -Wl,-Bsymbolic-functions -Wl,-Bsymbolic)
target_link_libraries(zerork_cfd_plugin_tester.x zerork_cfd_plugin)

if(ENABLE_OPENMP)
target_compile_definitions(zerork_cfd_plugin_tester.x PRIVATE USE_OMP)
target_link_libraries(zerork_cfd_plugin_tester.x OpenMP::OpenMP_CXX)
endif()

set(SPIFY_PLUGIN_TARGETS "zerork_cfd_plugin")
set(SPIFY_TESTER_TARGETS "zerork_cfd_plugin_tester.x")

set_target_properties(zerork_cfd_plugin PROPERTIES PUBLIC_HEADER  "zerork_cfd_plugin.h")
install(TARGETS zerork_cfd_plugin LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION include)
install(TARGETS zerork_cfd_plugin_tester.x)

add_spifyIFP_target(ZeroRKCFDPluginIFP ${SPIFY_PLUGIN_TARGETS})
add_spifyIFP_target(ZeroRKCFDPluginTesterIFP ${SPIFY_TESTER_TARGETS})

