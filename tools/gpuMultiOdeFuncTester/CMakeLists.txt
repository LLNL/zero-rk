
if(${SUNDIALS_VERSION} EQUAL "5")

add_executable(gpuMultiOdeFuncTester.x gpuMultiOdeFuncTester.cu)
add_executable(gpuMultiOdeFunc.x gpuMultiOdeFunc.cu)

target_link_libraries(gpuMultiOdeFunc.x zerork_cuda sundials_nvecserial sundials_nveccuda lapack blas)
target_link_libraries(gpuMultiOdeFuncTester.x zerork_cuda sundials_nvecserial sundials_nveccuda lapack blas)

endif()
