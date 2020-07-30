
include(FetchContent)
FetchContent_Declare(
  superlu_dist
  GIT_REPOSITORY https://github.com/xiaoyeli/superlu_dist
  GIT_TAG        v5.2.2 #TODO: can we move forward without fixing code
  GIT_SHALLOW    ON
)

FetchContent_GetProperties(superlu_dist)
if(NOT superlu_dist_POPULATED)
  FetchContent_Populate(superlu_dist)
  SET(enable_blaslib OFF CACHE INTERNAL "")
  SET(enable_tests OFF CACHE INTERNAL "")
  SET(enable_examples OFF CACHE INTERNAL "")
  SET(XSDK_ENABLE_Fortran OFF CACHE INTERNAL "")
  SET(TPL_PARMETIS_LIBRARIES ${parmetis_prefix}/lib/libparmetis.a CACHE INTERNAL "")
  SET(TPL_PARMETIS_INCLUDE_DIRS "${parmetis_prefix}/include" CACHE INTERNAL "")
  add_subdirectory(
    ${superlu_dist_SOURCE_DIR}
    ${superlu_dist_BINARY_DIR}
    EXCLUDE_FROM_ALL)
  #set_property(DIRECTORY ${superlu_dist_SOURCE_DIR} PROPERTY COMPILE_DEFINITIONS "METIS_EXPORT")
endif()

#TODO: why do we have to add fopenmp here (instead of depending on cmake to figure this out)
#TODO: likewise for libpthread
set_target_properties(superlu_dist PROPERTIES 
                      COMPILE_DEFINITIONS "METIS_EXPORT="
                      COMPILE_OPTIONS "-fopenmp"
                      INTERFACE_COMPILE_OPTIONS "-fopenmp"
                      INTERFACE_LINK_OPTIONS "-fopenmp"
                      INTERFACE_INCLUDE_DIRECTORIES ${superlu_dist_SOURCE_DIR}/SRC)
cmake_policy(SET CMP0079 NEW)
target_link_libraries(superlu_dist parmetis pthread)

