
#This module will make a header that contains
#the appropriate Windows definitions for
#creating/using a shared library
include(GenerateExportHeader)


set(COMMON_SRCS optionable.cpp
                zerork_flame_api.cpp
                zerork_flame_manager.cpp
                flame_params.cpp
                kinsol_functions.cpp
                sparse_matrix.cpp
                zerork_flame_api.cpp)
set(SPIFY_APPS  "zerork_flame_api_tester.x;zerork_flame_api")

add_library(zerork_flame_api SHARED ${COMMON_SRCS})
target_link_libraries(zerork_flame_api zerorkutilities mechanisminfo reactor
                      zerorktransport zerork superlu spify sundials_kinsol
                      sundials_nvecserial)
if(NOT WIN32 AND NOT APPLE)
target_link_options(zerork_flame_api PRIVATE -Wl,-Bsymbolic-functions -Wl,-Bsymbolic)
endif()
target_include_directories(zerork_flame_api PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

add_executable(zerork_flame_api_tester.x zerork_flame_api_tester.cpp)
target_link_libraries(zerork_flame_api_tester.x zerork_flame_api)

install(TARGETS zerork_flame_api PUBLIC_HEADER DESTINATION include
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION lib
        ARCHIVE DESTINATION lib)
install(TARGETS zerork_flame_api_tester.x
        RUNTIME DESTINATION bin)

if(ENABLE_MPI)
add_mpi_library(zerork_flame_api_mpi SHARED ${COMMON_SRCS} sparse_matrix_dist.cpp)
target_link_libraries(zerork_flame_api_mpi zerorkutilities mechanisminfo reactor
                      zerorktransport zerork superlu_dist superlu spify sundials_kinsol
                      sundials_nvecparallel)
if(NOT WIN32 AND NOT APPLE)
target_link_options(zerork_flame_api_mpi PRIVATE -Wl,-Bsymbolic-functions -Wl,-Bsymbolic)
endif()
target_include_directories(zerork_flame_api_mpi PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

add_mpi_executable(zerork_flame_api_tester_mpi.x zerork_flame_api_tester.cpp)
target_link_libraries(zerork_flame_api_tester_mpi.x zerork_flame_api_mpi)
set(SPIFY_APPS "${SPIFY_APPS};zerork_flame_api_tester_mpi.x;zerork_flame_api_mpi")

install(TARGETS zerork_flame_api_mpi PUBLIC_HEADER DESTINATION include
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION lib
        ARCHIVE DESTINATION lib)
install(TARGETS zerork_flame_api_tester_mpi.x
        RUNTIME DESTINATION bin)

endif()

add_spifyIFP_target(ZeroRKFlameAPITesterIFP ${SPIFY_APPS})

GENERATE_EXPORT_HEADER(zerork_flame_api
                       BASE_NAME zerork_flame
                       EXPORT_MACRO_NAME ZERORK_FLAME_EXPORTS
                       EXPORT_FILE_NAME zerork_flame_exports.h)

