
add_executable(randomStateGen.x randomStateGen.c)

target_link_libraries(randomStateGen.x m)

install(TARGETS randomStateGen.x)

