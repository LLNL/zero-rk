

add_library(zerork_equilibrium zerork_equilibrium_solver.cpp)

target_link_libraries(zerork_equilibrium zerork)
if(ZERORK_HAVE_CEQ)
  target_link_libraries(zerork_equilibrium zerork_ceq)
else()
  if(ZERORK_HAVE_CANTERA)
    target_link_libraries(zerork_equilibrium zerork_cantera)
  endif()
endif()

