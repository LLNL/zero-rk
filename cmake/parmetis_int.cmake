
include(FetchContent)
FetchContent_Declare(
  parmetis
  URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
  URL_HASH          MD5=f69c479586bf6bb7aff6a9bc0c739628
  PATCH_COMMAND patch -p0 < "${CMAKE_CURRENT_LIST_DIR}/parmetis-4.0.3.patch"
)

FetchContent_GetProperties(parmetis)
if(NOT parmetis_POPULATED)
  FetchContent_Populate(parmetis)

  add_subdirectory(
    ${parmetis_SOURCE_DIR}
    ${parmetis_BINARY_DIR})
endif()

