

include(FetchContent)
FetchContent_Declare(
  spify
  URL https://github.com/LLNL/spify/archive/7d43f1b375243344578057b80f0d2875c1bef295.tar.gz
)

FetchContent_GetProperties(spify)
if(NOT spify_POPULATED)
  FetchContent_Populate(spify)
  add_subdirectory(
    ${spify_SOURCE_DIR}
    ${spify_BINARY_DIR}
  )
endif()

