
#https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
#(accessed 20200513)

#TODO: Check if we have a clean source
#      and note in rev info
# https://stackoverflow.com/a/5139672

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
#TODO: Can we trigger this only on changes?
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
                  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                  OUTPUT_VARIABLE ZERORK_GIT_ID
                  ERROR_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${GIT_EXECUTABLE} show -s --format=%ci
                  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                  OUTPUT_VARIABLE ZERORK_GIT_COMMIT_TIMESTAMP
                  ERROR_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
                  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                  OUTPUT_VARIABLE ZERORK_GIT_BRANCH
                  ERROR_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

