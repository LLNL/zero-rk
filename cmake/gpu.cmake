
option(ENABLE_GPU "Enable CUDA libraries and applications." OFF)
cmake_dependent_option(ENABLE_MAGMA "Enable MAGMA GPU library" ON "ENABLE_GPU" OFF)
if(${ENABLE_GPU})
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14)
  set(CMAKE_CUDA_ARCHITECTURES "60;61;75")
  #N.B. this is for configuring external projects that use CMAKE_CUDA_ARCHITECTURES
  #     as semi-colon separated lists get munged when passing through ExternalProject_Add
  string(REPLACE ";" "|" CMAKE_CUDA_ARCHITECTURES_EXT "${CMAKE_CUDA_ARCHITECTURES}")

  set(ZERORK_HAVE_MAGMA OFF)
  if(${ENABLE_MAGMA})
    include(cmake/magma.cmake)
    set(ZERORK_HAVE_MAGMA ON)
  endif()
endif()

