

set(superlu_system_working ON)
if(EXISTS ${SYSTEM_SUPERLU_ROOT})
  set(superlu_prefix ${SYSTEM_SUPERLU_ROOT})
  find_library(superlu superlu
      PATHS ${superlu_prefix}
      PATH_SUFFIXES lib lib64)
  if(superlu-NOTFOUND)
    message(error "Couldn't find superlu_${LIBRARY}")
    set(superlu_system_working OFF)
  endif()
else()
  set(superlu_prefix ${CMAKE_CURRENT_BINARY_DIR}/external/superlu)
endif()

if((NOT EXISTS ${superlu_prefix}) OR (NOT ${superlu_system_working}))
  message(STATUS "Building: SuperLU...")
  configure_file(
	${CMAKE_CURRENT_LIST_DIR}/superlu.CMakeLists.txt
	${CMAKE_BINARY_DIR}/external/superlu-build/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/superlu-build)
  if(result)
    message(FATAL_ERROR "CMake step for SuperLU failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/superlu-build)
  if(result)
    message(FATAL_ERROR "Build step for SuperLU failed: ${result}")
  endif()
endif()

add_library(superlu STATIC IMPORTED GLOBAL)
set_target_properties(superlu PROPERTIES
  IMPORTED_LOCATION ${superlu_prefix}/lib64/${CMAKE_STATIC_LIBRARY_PREFIX}superlu${CMAKE_STATIC_LIBRARY_SUFFIX}
  INTERFACE_INCLUDE_DIRECTORIES ${superlu_prefix}/include)

target_link_libraries(superlu INTERFACE lapack blas)

