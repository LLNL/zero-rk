
option(ENABLE_GPU "Enable CUDA libraries and applications." OFF)
cmake_dependent_option(ENABLE_MAGMA "Enable MAGMA GPU library" ON "ENABLE_GPU" OFF)
if(${ENABLE_GPU})
  set(CMAKE_HIP_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  enable_language(HIP)
  if(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
    set(CMAKE_HIP_ARCHITECTURES gfx90a gfx942)
  endif()
  string(REPLACE ";" "|" CMAKE_HIP_ARCHITECTURES_EXT "${CMAKE_HIP_ARCHITECTURES}")
  message("CMAKE_HIP_ARCHITECTURES_EXT: ${CMAKE_HIP_ARCHITECTURES_EXT}")

  set(ZERORK_HAVE_MAGMA OFF)
  if(${ENABLE_MAGMA})
    include(cmake/magma.cmake)
    set(ZERORK_HAVE_MAGMA ON)
  endif()
endif()

