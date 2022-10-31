
add_library(reactor variable_volume_reactor.cpp const_pressure_reactor.cpp counterflow_reactor.cpp reactor_base.cpp)

target_include_directories(reactor PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                          $<INSTALL_INTERFACE:include>
                                   PRIVATE ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(reactor zerork)

set(public_headers const_pressure_reactor.h counterflow_reactor.h
reactor_base.h reactor_constants.h
variable_volume_reactor.h)

set_target_properties(reactor PROPERTIES
 PUBLIC_HEADER  "${public_headers}")
install(TARGETS reactor
 LIBRARY DESTINATION lib
 ARCHIVE DESTINATION lib
 PUBLIC_HEADER DESTINATION include/reactor)
