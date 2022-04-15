
add_executable(functionTester.x functionTester.cpp)

target_link_libraries(functionTester.x zerork ckconverter lapack blas)
if(NOT WIN32)
target_link_libraries(functionTester.x dl m)
endif()

install(TARGETS functionTester.x
        RUNTIME DESTINATION bin)

