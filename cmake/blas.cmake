

if(DEFINED ENV{MKLROOT})
set(ZERORK_HAVE_MKL ON CACHE BOOL "")
add_library(blas INTERFACE IMPORTED GLOBAL)
set_target_properties(blas PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES $ENV{MKLROOT}/include)
#Below pulled from MKL Link Advisor 2018.0.
target_link_libraries(blas INTERFACE
  -Wl,--start-group
  $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a
  $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a
  $ENV{MKLROOT}/lib/intel64/libmkl_core.a 
  -Wl,--end-group
  pthread m dl)
get_target_property(BLAS_LIBRARIES blas INTERFACE_LINK_LIBRARIES)
else()
set(ZERORK_HAVE_MKL OFF CACHE BOOL "")
find_package(BLAS REQUIRED)
add_library(blas INTERFACE IMPORTED GLOBAL)
target_link_libraries(blas INTERFACE ${BLAS_LIBRARIES})
endif()


