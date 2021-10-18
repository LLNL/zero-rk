
if(${ZERORK_HAVE_MKL})
  add_library(lapack ALIAS mkl)
else()
  find_package(LAPACK REQUIRED)
  add_library(lapack INTERFACE IMPORTED GLOBAL)
  target_link_libraries(lapack INTERFACE ${LAPACK_LIBRARIES} blas)
endif()

