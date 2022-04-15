
add_library(ckconverter CKParser.cpp CKReader.cpp ckr_utils.cpp Reaction.cpp thermoFunctions.cpp writelog.cpp)

target_include_directories(ckconverter PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                              $<INSTALL_INTERFACE:include>)

set(public_headers CKParser.h CKReader.h Constituent.h
Element.h Reaction.h RxnSpecies.h Species.h
atomicWeightDB.h ckr_defs.h ckr_utils.h thermoFunctions.h
writelog.h)
set_target_properties(ckconverter PROPERTIES
 PUBLIC_HEADER  "${public_headers}")
install(TARGETS ckconverter
 LIBRARY DESTINATION lib
 ARCHIVE DESTINATION lib
 PUBLIC_HEADER DESTINATION include/CKconverter)

