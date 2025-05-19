

if(ZERORK_HAVE_CEQ)
  set(CEQ_SRCS ceq_cp.f90 ceq_dgdt.f90 ceq_g.f90 ceq_h.f90
               ceq_h2T.f90 ceq_ind_col.f90 ceq_input.f90
               ceq_linprog.f90 ceq_lss.f90 ceq_maxmin.f90
               ceq_min_pert.f90 ceq_ming.f90 ceq_norm.f90
               ceq_perturb.f90 ceq_red_con.f90 ceq_reorder.f90
               ceq_solve.f90 ceq_state.f90 ceq_system.f90
               ceq_types.f90)
  list(TRANSFORM CEQ_SRCS PREPEND ${CMAKE_BINARY_DIR}/external/CEQ/)
  add_library(zerork_ceq ${CEQ_SRCS} eq_interface.f90 zerork_ceq_manager.cpp)
  target_link_libraries(zerork_ceq zerork lapack)
  target_compile_definitions(zerork_ceq PUBLIC ZERORK_HAVE_CEQ)
endif()

