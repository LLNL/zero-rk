

set(INPUT_FILES ./flamesolver.py.in ./counterflowflamesolver.py.in ./counterflowpremixedflamesolver.py.in)

foreach(INPUT_FILE ${INPUT_FILES})
string(REGEX REPLACE "\.in$" "" OUTPUT_FILE ${INPUT_FILE})
configure_file(${INPUT_FILE} ${OUTPUT_FILE} @ONLY)
get_filename_component(OUT_DIR ${INPUT_FILE} DIRECTORY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${OUTPUT_FILE} DESTINATION lib/python/site-packages/${OUT_DIR})
endforeach()

