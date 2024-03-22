
#WIP: Example of more detailed check
#https://github.com/QMCPACK/qmcpack/blob/develop/CMake/FindIBMMASS.cmake

include(CheckIncludeFile)
add_library(zerork_vectormath INTERFACE)
#IBM MASS
if(NOT DEFINED ZERORK_HAVE_IBM_MASS)
  if(DEFINED ENV{MASSROOT})
    CHECK_INCLUDE_FILE("massv.h" ZERORK_HAVE_IBM_MASS -I$ENV{MASSROOT}/include)
    message(STATUS "Found IBM MASS in $ENV{MASSROOT}")
  endif()
  if(${ZERORK_HAVE_IBM_MASS}) 
    set(ZERORK_HAVE_IBM_MASS ON CACHE BOOL "" FORCE)
    mark_as_advanced(ZERORK_HAVE_IBM_MASS)
    set(ZERORK_MASSROOT "$ENV{MASSROOT}" CACHE FILEPATH "" FORCE)
    mark_as_advanced(ZERORK_MASSROOT)
  else()
    set(ZERORK_HAVE_IBM_MASS OFF CACHE BOOL "" FORCE)
    mark_as_advanced(ZERORK_HAVE_IBM_MASS)
  endif()
endif()

if(NOT ZERORK_EXP_LIBC)
  if(${ZERORK_HAVE_IBM_MASS})
    add_library(IBM::MASSV INTERFACE IMPORTED GLOBAL)
    target_include_directories(IBM::MASSV INTERFACE "${ZERORK_MASSROOT}/include")
    target_link_libraries(IBM::MASSV INTERFACE "${ZERORK_MASSROOT}/lib/libmassvp9.a")
    target_link_libraries(zerork_vectormath INTERFACE IBM::MASSV)
  endif()

  if(${ZERORK_HAVE_MKL})
    target_link_libraries(zerork_vectormath INTERFACE mkl)
  endif()
endif()

if( (NOT ZERORK_HAVE_IBM_MASS AND NOT ZERORK_HAVE_MKL) OR ZERORK_EXP_LIBC)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "^x86_64$")
    #use FMATH
    target_link_libraries(zerork_vectormath INTERFACE zerork)
  else()
    if(NOT ZERORK_EXP_LIBC)
      set(ZERORK_EXP_LIBC ON CACHE BOOL "" FORCE)
      message(STATUS "No available fast exponential function.  Using libc exponential.")
    endif()
  endif()
endif()

