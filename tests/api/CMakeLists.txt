
set(SRCS elemental_composition_gtest.cpp physical_constants_gtest.cpp mechanism_gtest.cpp)

foreach(TEST_SRC ${SRCS})
string(REPLACE .cpp .x TEST ${TEST_SRC})
zerork_add_gtests(${TEST} SOURCES ${TEST_SRC} LINK_LIBRARIES zerork)
endforeach()

