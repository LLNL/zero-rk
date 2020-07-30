

include(FetchContent)
FetchContent_Declare(
  spify
  GIT_REPOSITORY https://github.com/LLNL/spify
  GIT_TAG        v1.0.4
  GIT_SHALLOW    ON
)

FetchContent_GetProperties(spify)
if(NOT spify_POPULATED)
  FetchContent_Populate(spify)
  add_subdirectory(
    ${spify_SOURCE_DIR}
    ${spify_BINARY_DIR}
  )
endif()

