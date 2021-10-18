
if(NOT DEFINED ZERORK_HAVE_MKL)
  if(DEFINED ENV{MKLROOT})
    set(ZERORK_HAVE_MKL ON CACHE BOOL "" FORCE)
    mark_as_advanced(ZERORK_HAVE_MKL)
    set(ZERORK_MKLROOT "$ENV{MKLROOT}" CACHE FILEPATH "" FORCE)
    mark_as_advanced(ZERORK_MKLROOT)
  else()
    set(ZERORK_HAVE_MKL OFF CACHE BOOL "" FORCE)
    mark_as_advanced(ZERORK_HAVE_MKL)
  endif()
endif()

if(${ZERORK_HAVE_MKL})
  add_library(mkl INTERFACE IMPORTED GLOBAL)
  #These removals are to work around an issue with mkl `module`s on LLNL's
  #clusters that set environment variables which cmake acts upon on first configure
  #but does not refresh on re-configuration.  This can mean that different paths
  #are used depending on environment state (independent of configuration variables)
  #which can cause the build to fail if the user forgets to `module load mkl` when
  #(re-)building the code.
  list(REMOVE_ITEM CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "${ZERORK_MKLROOT}/include")
  list(REMOVE_ITEM CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES "${ZERORK_MKLROOT}/include")
  list(REMOVE_ITEM CMAKE_C_IMPLICIT_LINK_DIRECTORIES "${ZERORK_MKLROOT}/lib/intel64")
  list(REMOVE_ITEM CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "${ZERORK_MKLROOT}/lib/intel64")
  target_include_directories(mkl INTERFACE "${ZERORK_MKLROOT}/include")
  #Below pulled from MKL Link Advisor 2018.0.
  target_link_libraries(mkl INTERFACE
    -Wl,--start-group
    "${ZERORK_MKLROOT}/lib/intel64/libmkl_intel_lp64.a"
    "${ZERORK_MKLROOT}/lib/intel64/libmkl_sequential.a"
    "${ZERORK_MKLROOT}/lib/intel64/libmkl_core.a"
    -Wl,--end-group
    pthread m dl)
   add_library(blas ALIAS mkl)
else()
  find_package(BLAS REQUIRED)
  add_library(blas INTERFACE IMPORTED GLOBAL)
  target_link_libraries(blas INTERFACE ${BLAS_LIBRARIES})
endif()

