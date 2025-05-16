
#User needs to supply a cantera built with compatible C++ and yaml-cpp
#The following build command worked with Cantera 3.0.1 on LC
#scons build boost_inc_dir=$BOOST_ROOT/include python_package=none system_yamlcpp=n prefix=/install/path
#Then set CANTERA_ROOT to /install/path

option(ENABLE_CANTERA "Enable interface to CANTERA. User must supply cantera." OFF)
set(ZERORK_HAVE_CANTERA OFF)

if(ENABLE_CANTERA)
  if(EXISTS ${CANTERA_ROOT})
    find_library(CANTERA_LIBRARY cantera PATHS ${CANTERA_ROOT}/lib)
    if(CANTERA_LIBRARY)
      set(CMAKE_CXX_STANDARD 17)
      set(ZERORK_HAVE_CANTERA ON)
      add_library(cantera STATIC IMPORTED GLOBAL)
      set_target_properties(cantera PROPERTIES
 	    IMPORTED_LOCATION ${CANTERA_LIBRARY}
	    INTERFACE_INCLUDE_DIRECTORIES ${CANTERA_ROOT}/include)
      target_link_libraries(cantera INTERFACE blas)
    else()
      message(STATUS "Can't find Cantera library.  Disabling Cantera interface.")
    endif()
  else()
    message(STATUS "Can't find Cantera directory.  Disabling Cantera interface.")
  endif()
endif()

