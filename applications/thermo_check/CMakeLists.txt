

set(COMMON_SRCS thermo_parser.cpp thermo_fix.cpp
                janaf_thermo.cpp polynomial.c
                lapack_wrapper.c ThermoCheckIFP.cpp)

add_executable(thermo_check.x thermo_check.cpp ${COMMON_SRCS})

add_executable(thermo_check_only.x thermo_check_only.cpp 
               ${COMMON_SRCS})

target_link_libraries(thermo_check.x zerork spify zerorkutilities lapack)
target_link_libraries(thermo_check_only.x zerork spify zerorkutilities lapack)
add_spifyIFP_target(ThermoCheckIFP thermo_check.x thermo_check_only.x)

install(TARGETS thermo_check.x thermo_check_only.x
        RUNTIME DESTINATION bin)

