

add_executable(gpuMultiOdeFuncTester.x gpuMultiOdeFuncTester.hip)
add_executable(gpuMultiOdeFunc.x gpuMultiOdeFunc.hip)

target_link_libraries(gpuMultiOdeFunc.x zerork_cuda sundials_nvecserial sundials_nvechip lapack blas)
target_link_libraries(gpuMultiOdeFuncTester.x zerork_cuda sundials_nvecserial sundials_nvechip lapack blas)

