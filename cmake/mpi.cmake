

find_package(MPI REQUIRED COMPONENTS C CXX)

function(add_mpi_executable)
set(TARGET ${ARGV0})
add_executable(${ARGV}) #ARGV is all the arguments

target_include_directories(${TARGET} PRIVATE ${MPI_INCLUDE_PATH})

target_link_libraries(${TARGET} ${MPI_LIBRARIES})

set(FLAGS "-DZERORK_MPI ${MPI_COMPILE_FLAGS}")
set_target_properties(${TARGET} PROPERTIES
    COMPILE_FLAGS "${FLAGS}")

if(MPI_LINK_FLAGS)
  set_target_properties(${TARGET} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()#
endfunction()

function(add_mpi_library)
set(TARGET ${ARGV0})
add_library(${ARGV}) #ARGV is all the arguments

target_include_directories(${TARGET} PRIVATE ${MPI_INCLUDE_PATH})

target_link_libraries(${TARGET} ${MPI_LIBRARIES})

set(FLAGS "-DZERORK_MPI ${MPI_COMPILE_FLAGS}")
set_target_properties(${TARGET} PROPERTIES
    COMPILE_FLAGS "${FLAGS}")

if(MPI_LINK_FLAGS)
  set_target_properties(${TARGET} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()#
endfunction()

