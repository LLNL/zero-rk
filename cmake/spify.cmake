

include(FetchContent)
FetchContent_Declare(
  spify
  URL https://github.com/LLNL/spify/archive/refs/tags/v1.0.11.tar.gz
)

FetchContent_GetProperties(spify)
if(NOT spify_POPULATED)
  FetchContent_Populate(spify)
  add_subdirectory(
    ${spify_SOURCE_DIR}
    ${spify_BINARY_DIR}
  )
endif()

