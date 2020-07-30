
set(SRCS test01_const_pressure_derivative.cpp 
  test02_const_pressure_equilibrium.cpp
  test03_const_pressure_jacobian.cpp
  test05_const_pressure_time_history.cpp
  test06_const_pressure_derivative.cpp)

foreach(TEST_SRC ${SRCS})
string(REPLACE .cpp .x TEST ${TEST_SRC})
zerork_add_test(${TEST}
    SOURCES ${TEST_SRC}
    LINK_LIBRARIES reactor zerorktransport zerork
                   sundials_cvode sundials_nvecserial
    ARGUMENTS 10 -1)
endforeach()

#Multiple source files for test04
zerork_add_test(test04_const_pressure_jacobian.x
    SOURCES test04_const_pressure_jacobian.cpp test04_utilities.cpp
    LINK_LIBRARIES reactor zerorktransport zerork
                   sundials_cvode sundials_nvecserial
    ARGUMENTS 10 -1)

