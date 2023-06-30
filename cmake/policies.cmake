
# enable MSVC_RUNTIME_LIBRARY target property
# see https://cmake.org/cmake/help/latest/policy/CMP0091.html
if(POLICY CMP0091)
  cmake_policy(SET CMP0091 NEW)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

# avoid warnings about DOWNLOAD_EXTRACT_TIMESTAMP in ExternalProject_Add
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

