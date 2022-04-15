
add_executable(randomStateGen.x randomStateGen.cpp)

target_link_libraries(randomStateGen.x zerorkutilities)
if(NOT WIN32)
target_link_libraries(randomStateGen.x m)
endif()

install(TARGETS randomStateGen.x
        RUNTIME DESTINATION bin)

