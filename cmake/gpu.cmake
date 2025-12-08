
option(ENABLE_GPU "Enable CUDA/HIP libraries and applications." OFF)
cmake_dependent_option(ENABLE_MAGMA "Enable MAGMA GPU library" ON "ENABLE_GPU" OFF)

if(ENABLE_GPU)
  set(CMAKE_HIP_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  enable_language(HIP)

  #Determine arch from env variables
  if(DEFINED ENV{MI300A})
    message(STATUS "env var MI300A set, using MI300A arch with gfx942")
    set(HIP_ARCHITECTURES gfx942)
    set(CMAKE_HIP_ARCHITECTURES gfx942)

  elseif(DEFINED ENV{MI250})
    message(STATUS "env var MI250 set, using MI250 arch with gfx90a")
    set(HIP_ARCHITECTURES gfx90a)
    set(CMAKE_HIP_ARCHITECTURES gfx90a)

  #If arch is not set from env variable:
  elseif(DEFINED ENV{HIPCC_COMPILE_FLAGS_APPEND})
    if("$ENV{HIPCC_COMPILE_FLAGS_APPEND}" MATCHES "gfx942")
      message(STATUS "using gfx942 from HIPCC_COMPILE_FLAGS_APPEND")
      set(HIP_ARCHITECTURES gfx942)
      set(CMAKE_HIP_ARCHITECTURES gfx942)

    elseif("$ENV{HIPCC_COMPILE_FLAGS_APPEND}" MATCHES "gfx90a")
      message(STATUS "using gfx90a from HIPCC_COMPILE_FLAGS_APPEND")
      set(HIP_ARCHITECTURES gfx90a)
      set(CMAKE_HIP_ARCHITECTURES gfx90a)

    else()
      message(STATUS "could not detect gfx90a or gfx942 architecture from hipcc compile flag, defaulting to gfx942")
      set(HIP_ARCHITECTURES gfx942)
      set(CMAKE_HIP_ARCHITECTURES gfx942)
    endif()

  else()
    message(STATUS "no env variables were given, CMAKE defaults to the MI300A")
    set(HIP_ARCHITECTURES gfx942)
    set(CMAKE_HIP_ARCHITECTURES gfx942)
  endif()


  string(REPLACE ";" "|" CMAKE_HIP_ARCHITECTURES_EXT "${CMAKE_HIP_ARCHITECTURES}")

  set(ZERORK_HAVE_MAGMA OFF)
  if(ENABLE_MAGMA)
    include(cmake/magma.cmake)
    set(ZERORK_HAVE_MAGMA ON)
  endif()

endif()
