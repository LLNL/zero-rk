
set(SRCS big_molecule_gtest.cpp non_integer_gtest.cpp
   plog_gtest.cpp sri_gtest.cpp troe_gtest.cpp rord_propane_gtest.cpp
   ford_propane_gtest.cpp)

foreach(TEST_SRC ${SRCS})
string(REPLACE .cpp .x TEST ${TEST_SRC})
zerork_add_gtests(${TEST} SOURCES ${TEST_SRC} LINK_LIBRARIES zerork)
endforeach()

zerork_add_test(plog_test.x
    SOURCES plog_test.cpp
    LINK_LIBRARIES zerork zerorkutilities
    ARGUMENTS 1000 -1 1.0e-20 1.0e-8)
#add_test(NAME plog_test.x COMMAND plog_test.x 1000 -1 1.0e-20 1.0e-8)

   
