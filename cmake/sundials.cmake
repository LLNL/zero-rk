
set(SUNDIALS_VERSION "5" CACHE STRING "Sundials major version interface")
set_property(CACHE SUNDIALS_VERSION PROPERTY STRINGS "2;3;4;5")

cmake_dependent_option(ZERORK_ENABLE_SUNDIALS_LAPACK 
                       "Enable LAPACK in Sundials; only for Sundials-5 on Linux/Mac"
                       ON "SUNDIALS_VERSION EQUAL 5;NOT WIN32" OFF)

set(sundials_libs nvecserial nvecparallel cvode cvodes kinsol ida arkode)

if(${ZERORK_ENABLE_SUNDIALS_LAPACK})
  set(sundials_libs "${sundials_libs};sunlinsollapackdense")
endif()

set(SUNDIALS_IFACE_NUM ${SUNDIALS_VERSION})
if(${SUNDIALS_VERSION} EQUAL "5")
set(SUNDIALS_IFACE_NUM "4")
if(${ENABLE_GPU})
set(sundials_libs "${sundials_libs};nveccuda")
set(MAGMA_ARGS "")
if(${ZERORK_HAVE_MAGMA})
set(sundials_libs "${sundials_libs};sunlinsolmagmadense;sunmatrixmagmadense")
set(MAGMA_ARGS "-DENABLE_MAGMA=ON -DMAGMA_DIR=${magma_prefix}")
endif() #MAGMA
endif() #GPU
endif() #Sundials5

set(sundials_system_working ON)
if(EXISTS ${SYSTEM_SUNDIALS_ROOT})
  set(sundials_prefix ${SYSTEM_SUNDIALS_ROOT})
  foreach(LIBRARY ${sundials_libs})
    find_library(sundials_${LIBRARY}
      ${LIBRARY}
      PATHS ${sundials_prefix}
      PATH_SUFFIXES lib lib64)
    if(sundials_${LIBRARY}-NOTFOUND)
      message(error "Couldn't find sundials_${LIBRARY}")
      set(sundials_system_working OFF)
    endif()
  endforeach()
else()
  set(sundials_prefix ${CMAKE_CURRENT_BINARY_DIR}/external/sundials${SUNDIALS_VERSION}/${ZERORK_EXTERNALS_BUILD_TYPE})
endif()

if((NOT EXISTS ${sundials_prefix}) OR (NOT ${sundials_system_working}))
  message(STATUS "Building: Sundials-${SUNDIALS_VERSION}...")
  configure_file(
	${CMAKE_CURRENT_LIST_DIR}/sundials${SUNDIALS_VERSION}.CMakeLists.txt
	${CMAKE_BINARY_DIR}/external/sundials${SUNDIALS_VERSION}-build/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/sundials${SUNDIALS_VERSION}-build)
  if(result)
    message(FATAL_ERROR "CMake step for Sundials-${SUNDIALS_VERSION} failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build . --config ${ZERORK_EXTERNALS_BUILD_TYPE}
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/sundials${SUNDIALS_VERSION}-build)
  if(result)
    message(FATAL_ERROR "Build step for Sundials-${SUNDIALS_VERSION} failed: ${result}")
  endif()
endif()

set(sundials_lib_dir ${sundials_prefix}/lib64)
if(NOT EXISTS ${sundials_lib_dir})
  set(sundials_lib_dir ${sundials_prefix}/lib)
  if(NOT EXISTS ${sundials_lib_dir})
    message(FATAL_ERROR "Couldn't find sundials library directory.")
  endif()
endif()

foreach(LIBRARY ${sundials_libs})
  add_library(sundials_${LIBRARY} STATIC IMPORTED GLOBAL)
  set_target_properties(sundials_${LIBRARY} PROPERTIES
  IMPORTED_LOCATION ${sundials_lib_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}sundials_${LIBRARY}${CMAKE_STATIC_LIBRARY_SUFFIX}
  INTERFACE_INCLUDE_DIRECTORIES ${sundials_prefix}/include
  INTERFACE_COMPILE_DEFINITIONS SUNDIALS${SUNDIALS_IFACE_NUM})
  if(${ZERORK_HAVE_MAGMA})
    target_link_libraries(sundials_${LIBRARY} INTERFACE magma)
  endif()
  if(${ZERORK_ENABLE_SUNDIALS_LAPACK})
    set_property(TARGET sundials_${LIBRARY} APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS ZERORK_HAVE_SUNDIALS_LAPACK)
  endif()
  ##target_compile_definitions(sundials_${LIBRARY} INTERFACE SUNDIALS${SUNDIALS_IFACE_NUM} JACOBIAN_SIZE_LONG_INT)
endforeach()
