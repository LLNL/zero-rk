
option(ENABLE_GPU "Enable CUDA libraries and applications." OFF)
cmake_dependent_option(ENABLE_MAGMA "Enable MAGMA GPU library" ON "ENABLE_GPU" OFF)
if(${ENABLE_GPU})
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14)
  set(CMAKE_CUDA_ARCHITECTURES 60 61 75)

  set(ZERORK_HAVE_MAGMA OFF)
  if(${ENABLE_MAGMA})
    include(cmake/magma.cmake)
    set(ZERORK_HAVE_MAGMA ON)
  endif()
endif()

