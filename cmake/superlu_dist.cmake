

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
  set(superlu_dist_prefix ${CMAKE_CURRENT_BINARY_DIR}/external/superlu_dist/${ZERORK_EXTERNALS_BUILD_TYPE})
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
  if(WIN32)
    #superlu_dist/SRC/util.c includes unistd.h for no reason.  This is so we don't have to patch the code
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/external/superlu_dist-build/superlu_dist-prefix/src/superlu_dist-build/SRC)
    file(TOUCH ${CMAKE_BINARY_DIR}/external/superlu_dist-build/superlu_dist-prefix/src/superlu_dist-build/SRC/unistd.h)
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build . --config ${ZERORK_EXTERNALS_BUILD_TYPE}
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/superlu_dist-build)
  if(result)
    message(FATAL_ERROR "Build step for SuperLU_DIST failed: ${result}")
  endif()
endif()

set(superlu_dist_lib_dir ${superlu_dist_prefix}/lib64)
if(NOT EXISTS ${superlu_dist_lib_dir})
  set(superlu_dist_lib_dir ${superlu_dist_prefix}/lib)
  if(NOT EXISTS ${superlu_dist_lib_dir})
    message(FATAL_ERROR "Couldn't find superlu_dist library directory.")
  endif()
endif()

add_library(superlu_dist STATIC IMPORTED GLOBAL)
set_target_properties(superlu_dist PROPERTIES
  IMPORTED_LOCATION ${superlu_dist_lib_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}superlu_dist${CMAKE_STATIC_LIBRARY_SUFFIX}
  INTERFACE_INCLUDE_DIRECTORIES ${superlu_dist_prefix}/include)

target_link_libraries(superlu_dist INTERFACE lapack blas)

