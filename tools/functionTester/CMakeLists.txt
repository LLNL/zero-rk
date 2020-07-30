
add_executable(functionTester.x functionTester.cpp)

target_link_libraries(functionTester.x zerork ckconverter lapack blas dl m)

install(TARGETS functionTester.x)

