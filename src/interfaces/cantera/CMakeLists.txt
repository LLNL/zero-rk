

if(ZERORK_HAVE_CANTERA)
  add_library(zerork_cantera zerork_cantera_manager.cpp)
  target_link_libraries(zerork_cantera cantera zerork)
  target_compile_definitions(zerork_cantera PUBLIC ZERORK_HAVE_CANTERA)
endif()

