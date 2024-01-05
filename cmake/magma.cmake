

set(magma_system_working ON)
if(EXISTS ${SYSTEM_MAGMA_ROOT})
  set(magma_prefix ${SYSTEM_MAGMA_ROOT})
  find_library(magma magma
      PATHS ${magma_prefix}
      PATH_SUFFIXES lib lib64)
  if(magma-NOTFOUND)
    message(error "Couldn't find magma_${LIBRARY}")
    set(magma_system_working OFF)
  endif()
else()
  set(magma_prefix ${CMAKE_CURRENT_BINARY_DIR}/external/magma/${ZERORK_EXTERNALS_BUILD_TYPE})
endif()

if((NOT EXISTS ${magma_prefix}) OR (NOT ${magma_system_working}))
  message(STATUS "Building: Magma...")
  set(GPU_TARGET_LIST "${CMAKE_CUDA_ARCHITECTURES}")
  list(TRANSFORM GPU_TARGET_LIST PREPEND sm_)
  string(REPLACE ";" " " GPU_TARGET "${GPU_TARGET_LIST}")
  configure_file(
	${CMAKE_CURRENT_LIST_DIR}/magma.CMakeLists.txt
	${CMAKE_BINARY_DIR}/external/magma-build/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/magma-build)
  if(result)
    message(FATAL_ERROR "CMake step for Magma failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build . --config ${ZERORK_EXTERNALS_BUILD_TYPE}
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/magma-build)
  if(result)
    message(FATAL_ERROR "Build step for Magma failed: ${result}")
  endif()
endif()

set(magma_lib_dir ${magma_prefix}/lib64)
if(NOT EXISTS ${magma_lib_dir})
  set(magma_lib_dir ${magma_prefix}/lib)
  if(NOT EXISTS ${magma_lib_dir})
    message(FATAL_ERROR "Couldn't find magma library directory.")
  endif()
endif()

add_library(magma STATIC IMPORTED GLOBAL)
set_target_properties(magma PROPERTIES
  IMPORTED_LOCATION ${magma_lib_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}magma${CMAKE_STATIC_LIBRARY_SUFFIX}
  INTERFACE_INCLUDE_DIRECTORIES "${magma_prefix}/include;${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
  INTERFACE_LINK_DIRECTORIES "${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES}"
  INTERFACE_COMPILE_DEFINITIONS ZERORK_HAVE_MAGMA)


find_package(OpenMP)
target_link_libraries(magma INTERFACE cusparse lapack blas OpenMP::OpenMP_CXX)

