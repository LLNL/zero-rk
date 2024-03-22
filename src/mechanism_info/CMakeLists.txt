
add_library(mechanisminfo mechanism_info.cpp)

target_include_directories(mechanisminfo PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                                  $<INSTALL_INTERFACE:include>
                                   PRIVATE ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(mechanisminfo PUBLIC zerork)

set(public_headers mechanism_info.h)

set_target_properties(mechanisminfo PROPERTIES
 PUBLIC_HEADER  "${public_headers}")
install(TARGETS mechanisminfo
 LIBRARY DESTINATION lib
 ARCHIVE DESTINATION lib
 PUBLIC_HEADER DESTINATION include/mechanism_info)

