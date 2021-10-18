
#N.B. CMAKE_CURRENT_BINARY_DIR needed to include these files
configure_file(zerork_version_info.h.in zerork_version_info.h)
configure_file(zerork_conf.h.in zerork_conf.h)

add_library(zerork element.cpp species.cpp mechanism.cpp utilities.cpp
            nasa_poly.cpp info_net.cpp rate_const.cpp perf_net.cpp
            fast_exps.cpp external_funcs_blank.cpp plog_reaction.cpp
            non_integer_reaction_network.cpp constants_api.cpp 
            elemental_composition.cpp impls/elemental_composition_impl.cpp)

#TODO: Cleanup include directories (e.g. "zerork/rate_const.h" vs. "rate_const.h")
target_include_directories(zerork PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../>
                                         $<INSTALL_INTERFACE:include>)

set(public_headers atomicMassDB.h constants.h constants_api.h
   element.h elemental_composition.h external_funcs.h fast_exps.h
   info_net.h mechanism.h nasa_poly.h non_integer_reaction_network.h
   perf_net.h plog_reaction.h rate_const.h simple_network.h
   species.h utilities.h)
target_link_libraries(zerork PUBLIC ckconverter dl m)
if(${ZERORK_HAVE_MKL})
target_link_libraries(zerork PUBLIC mkl) #for vmdExp in fast_exp.cpp
endif()

set_target_properties(zerork PROPERTIES
 PUBLIC_HEADER  "${public_headers}")
install(TARGETS zerork
 LIBRARY DESTINATION lib
 PUBLIC_HEADER DESTINATION include/zerork)

