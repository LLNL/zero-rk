

set(INPUT_FILES setup.py.in
                src/rate_opt/config.py.in)

foreach(INPUT_FILE ${INPUT_FILES})
string(REGEX REPLACE "\.in$" "" OUTPUT_FILE ${INPUT_FILE})
configure_file(${INPUT_FILE} ${OUTPUT_FILE} @ONLY)
get_filename_component(OUT_DIR ${INPUT_FILE} DIRECTORY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${OUTPUT_FILE} DESTINATION share/zerork/python/${OUT_DIR})
endforeach()

install(DIRECTORY src DESTINATION share/zerork/python
        REGEX "\.in$" EXCLUDE)

install(FILES pyproject.toml DESTINATION share/zerork/python)

