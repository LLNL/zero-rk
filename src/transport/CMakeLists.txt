
add_library(zerorktransport binary_collision.cpp collision_integrals.cpp
                            mass_transport_factory.cpp constant_lewis.cpp
                            mix_avg.cpp mix_avg_soret.cpp flexible_transport.cpp)

target_include_directories(zerorktransport PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                                  $<INSTALL_INTERFACE:include>
                                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(zerorktransport zerorkutilities)

set(public_headers binary_collision.h collision_integrals.h
constant_lewis.h mass_transport_factory.h mix_avg.h mix_avg_soret.h flexible_transport.h)

set_target_properties(zerorktransport PROPERTIES
 PUBLIC_HEADER  "${public_headers}")
install(TARGETS zerorktransport
 LIBRARY DESTINATION lib
 ARCHIVE DESTINATION lib
 PUBLIC_HEADER DESTINATION include/transport)
