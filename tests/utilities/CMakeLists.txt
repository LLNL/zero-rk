

set(SRCS file_utilities_gtest.cpp math_utilities_gtest.cpp
         string_utilities_gtest.cpp)

foreach(TEST_SRC ${SRCS})
string(REPLACE .cpp .x TEST ${TEST_SRC})
zerork_add_gtests(${TEST} SOURCES ${TEST_SRC} LINK_LIBRARIES zerorkutilities)
endforeach()

#Only copy data once
add_dependencies(file_utilities_gtest.x math_utilities_gtest.x)
add_custom_command(TARGET math_utilities_gtest.x  PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_CURRENT_SOURCE_DIR}/data
                   $<TARGET_FILE_DIR:math_utilities_gtest.x>/data)

