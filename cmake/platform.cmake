
include(CheckSymbolExists)
check_symbol_exists(aligned_alloc "stdlib.h" HAVE_ALIGNED_ALLOC)

if(WIN32)
set(DEVNUL "NUL")
else()
set(DEVNUL "/dev/null")
endif()

