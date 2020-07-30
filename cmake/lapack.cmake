
if(DEFINED ENV{MKLROOT})
add_library(lapack INTERFACE IMPORTED GLOBAL)
set_target_properties(lapack PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES $ENV{MKLROOT}/include)
#Below pulled from MKL Link Advisor 2018.0.
target_link_libraries(lapack INTERFACE
  -Wl,--start-group
  $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a
  $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a
  $ENV{MKLROOT}/lib/intel64/libmkl_core.a 
  -Wl,--end-group
  blas
  pthread m dl)
get_target_property(LAPACK_LIBRARIES lapack INTERFACE_LINK_LIBRARIES)
else()
find_package(LAPACK REQUIRED)
add_library(lapack INTERFACE IMPORTED GLOBAL)
target_link_libraries(lapack INTERFACE ${LAPACK_LIBRARIES} blas)
endif()

