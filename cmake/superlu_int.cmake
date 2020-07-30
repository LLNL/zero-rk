
include(FetchContent)
FetchContent_Declare(
  superlu
  GIT_REPOSITORY https://github.com/xiaoyeli/superlu
  #GIT_TAG        v5.2.1
  GIT_SHALLOW    ON
)

FetchContent_GetProperties(superlu)
if(NOT superlu_POPULATED)
  FetchContent_Populate(superlu)
  SET(enable_blaslib OFF CACHE INTERNAL "")
  SET(XSDK_ENABLE_Fortran OFF CACHE INTERNAL "")

  add_subdirectory(
    ${superlu_SOURCE_DIR}
    ${superlu_BINARY_DIR}
    )
endif()

#TODO: Consider ExternalProject if trouble occurs with mixing cmake projects
#include(ExternalProject)
#ExternalProject_Add(
#  superlu_pkg
#  GIT_REPOSITORY https://github.com/xiaoyeli/superlu
#  GIT_TAG        v5.2.1
#  PREFIX ${CMAKE_BINARY_DIR}/external
#  CMAKE_ARGS -DCMAKE_enable_blaslib=OFF -DCMAKE_INSTALL_PREFIX="."
#)

#add_library(timestwo STATIC IMPORTED GLOBAL)
#set_target_properties(timestwo PROPERTIES
#  IMPORTED_LOCATION ${BINARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}timestwo${CMAKE_STATIC_LIBRARY_SUFFIX}
#  INTERFACE_INCLUDE_DIRECTORIES ${BINARY_DIR}/include)

#option(enable_blaslib   "Build the CBLAS library" ${enable_blaslib_DEFAULT})
#option(enable_matlabmex "Build the Matlab mex library" OFF)
#option(enable_tests     "Build tests" ON)
#option(enable_doc       "Build doxygen documentation" OFF)
#option(enable_single    "Enable single precision library" ON)
#option(enable_double    "Enable double precision library" ON)
#option(enable_complex   "Enable complex precision library" ON)
#option(enable_complex16 "Enable complex16 precision library" ON)

