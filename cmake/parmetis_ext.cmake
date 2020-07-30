
set(parmetis_system_working ON)
if(EXISTS ${SYSTEM_PARMETIS_ROOT})
  set(parmetis_prefix ${SYSTEM_PARMETIS_ROOT})
  find_library(parmetis parmetis
      PATHS ${parmetis_prefix}
      PATH_SUFFIXES lib lib64)
  if(parmetis-NOTFOUND)
    message(error "Couldn't find parmetis_${LIBRARY}")
    set(parmetis_system_working OFF)
  endif()
else()
  set(parmetis_prefix ${CMAKE_CURRENT_BINARY_DIR}/external/parmetis)
endif()

if((NOT EXISTS ${parmetis_prefix}) OR (NOT ${parmetis_system_working}))
  message(STATUS "Building: parmetis ...")
  configure_file(
	${CMAKE_CURRENT_LIST_DIR}/parmetis.CMakeLists.txt
	${CMAKE_BINARY_DIR}/external/parmetis-build/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/parmetis-build)
  if(result)
    message(FATAL_ERROR "CMake step for parmetis failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/parmetis-build)
  if(result)
    message(FATAL_ERROR "Build step for parmetis failed: ${result}")
  endif()
endif()

add_library(parmetis STATIC IMPORTED GLOBAL)
set_target_properties(parmetis PROPERTIES
  IMPORTED_LOCATION ${parmetis_prefix}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}parmetis${CMAKE_STATIC_LIBRARY_SUFFIX}
  INTERFACE_INCLUDE_DIRECTORIES ${parmetis_prefix}/include)

add_library(metis STATIC IMPORTED GLOBAL)
set_target_properties(metis PROPERTIES
  IMPORTED_LOCATION ${parmetis_prefix}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}metis${CMAKE_STATIC_LIBRARY_SUFFIX}
  INTERFACE_INCLUDE_DIRECTORIES ${parmetis_prefix}/include)

target_link_libraries(parmetis INTERFACE metis)
