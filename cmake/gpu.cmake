
option(ENABLE_GPU "Enable CUDA libraries and applications." OFF)
cmake_dependent_option(ENABLE_MAGMA "Enable MAGMA GPU library" ON "ENABLE_GPU" OFF)
if(${ENABLE_GPU})
  set(CMAKE_HIP_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  enable_language(HIP)
  #set(CMAKE_CUDA_STANDARD 14)
  #set(CMAKE_CUDA_ARCHITECTURES 60 61 75)
  set(HIP_ARCHITECTURES gfx90a)
  set(CMAKE_HIP_ARCHITECTURES gfx90a)
  string(REPLACE ";" "|" CMAKE_HIP_ARCHITECTURES_EXT "${CMAKE_HIP_ARCHITECTURES}")

  set(ZERORK_HAVE_MAGMA OFF)
  if(${ENABLE_MAGMA})
    include(cmake/magma.cmake)
    set(ZERORK_HAVE_MAGMA ON)
  endif()
endif()

