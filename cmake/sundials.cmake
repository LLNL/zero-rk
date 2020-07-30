
set(SUNDIALS_IFACE_NUM ${SUNDIALS_VERSION})
if(${SUNDIALS_VERSION} EQUAL "5")
set(SUNDIALS_IFACE_NUM "4")
endif()

set(sundials_libs nvecserial nvecparallel cvode cvodes kinsol arkode)

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
  set(sundials_prefix ${CMAKE_CURRENT_BINARY_DIR}/external/sundials${SUNDIALS_VERSION})
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
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/sundials${SUNDIALS_VERSION}-build)
  if(result)
    message(FATAL_ERROR "Build step for Sundials-${SUNDIALS_VERSION} failed: ${result}")
  endif()
endif()

foreach(LIBRARY ${sundials_libs})
  add_library(sundials_${LIBRARY} STATIC IMPORTED GLOBAL)
  set_target_properties(sundials_${LIBRARY} PROPERTIES
  IMPORTED_LOCATION ${sundials_prefix}/lib64/${CMAKE_STATIC_LIBRARY_PREFIX}sundials_${LIBRARY}${CMAKE_STATIC_LIBRARY_SUFFIX}
  INTERFACE_INCLUDE_DIRECTORIES ${sundials_prefix}/include
  INTERFACE_COMPILE_DEFINITIONS SUNDIALS${SUNDIALS_IFACE_NUM})
endforeach()


