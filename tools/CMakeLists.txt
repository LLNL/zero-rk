
add_subdirectory(functionTester)
add_subdirectory(randomStateGen)

if(ENABLE_GPU)
add_subdirectory(gpuMultiOdeFuncTester)
endif()

