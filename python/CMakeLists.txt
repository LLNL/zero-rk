

find_package (Python3 COMPONENTS Interpreter)

if(Python3_FOUND)
add_subdirectory(flame_solvers)
add_subdirectory(rate_optimization)
endif()

