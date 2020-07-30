

set(superlu_dist_system_working ON)
if(EXISTS ${SYSTEM_SUPERLU_DIST_ROOT})
  set(superlu_dist_prefix ${SYSTEM_SUPERLU_DIST_ROOT})
  find_library(superlu_dist superlu_dist
      PATHS ${superlu_dist_prefix}
      PATH_SUFFIXES lib lib64)
  if(superlu_dist-NOTFOUND)
    message(error "Couldn't find superlu_dist_${LIBRARY}")
    set(superlu_dist_system_working OFF)
  endif()
else()
  set(superlu_dist_prefix ${CMAKE_CURRENT_BINARY_DIR}/external/superlu_dist)
endif()

if((NOT EXISTS ${superlu_dist_prefix}) OR (NOT ${superlu_dist_system_working}))
  message(STATUS "Building: SuperLU_DIST...")
  configure_file(
	${CMAKE_CURRENT_LIST_DIR}/superlu_dist.CMakeLists.txt
	${CMAKE_BINARY_DIR}/external/superlu_dist-build/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/superlu_dist-build)
  if(result)
    message(FATAL_ERROR "CMake step for SuperLU_DIST failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/superlu_dist-build)
  if(result)
    message(FATAL_ERROR "Build step for SuperLU_DIST failed: ${result}")
  endif()
endif()

add_library(superlu_dist STATIC IMPORTED GLOBAL)
set_target_properties(superlu_dist PROPERTIES
  IMPORTED_LOCATION ${superlu_dist_prefix}/lib64/${CMAKE_STATIC_LIBRARY_PREFIX}superlu_dist${CMAKE_STATIC_LIBRARY_SUFFIX}
  INTERFACE_INCLUDE_DIRECTORIES ${superlu_dist_prefix}/include)

target_link_libraries(superlu_dist INTERFACE parmetis lapack blas pthread)

#TODO: This is clunky and possibly brittle
find_package(OpenMP)
if(OPENMP_FOUND)
  set_target_properties(superlu_dist PROPERTIES 
                        INTERFACE_COMPILE_OPTIONS "${OpenMP_C_FLAGS}"
                        INTERFACE_LINK_OPTIONS "${OpenMP_C_FLAGS}")
endif()

#TODO: why do we have to add fopenmp here (instead of depending on cmake to figure this out)
#TODO: likewise for libpthread
#set_target_properties(superlu_dist PROPERTIES 
#                      COMPILE_DEFINITIONS "METIS_EXPORT="
#                      COMPILE_OPTIONS "-fopenmp"
#                      INTERFACE_COMPILE_OPTIONS "-fopenmp"
#                      INTERFACE_LINK_OPTIONS "-fopenmp"
#                      INTERFACE_INCLUDE_DIRECTORIES ${superlu_dist_SOURCE_DIR}/SRC)

