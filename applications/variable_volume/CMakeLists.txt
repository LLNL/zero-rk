
set(TARGET variable_volume.x)
set(SRCS VariableVolumeIFP.cpp
         variable_volume.cpp complete_solver.cpp
         volume_functions.cpp sparse_matrix.cpp
         user_functions.cpp)

add_executable(${TARGET} ${SRCS})
target_link_libraries(${TARGET} mechanisminfo reactor zerork zerorkutilities
                      sundials_nvecserial sundials_cvode superlu spify)

install(TARGETS variable_volume.x
        RUNTIME DESTINATION bin)
add_spifyIFP_target(VariableVolumeIFP ${TARGET})

