set(AMICI_DIR "/home/dweindl/src/AMICI.github" CACHE PATH "AMICI root directory")

include_directories(
    src
    ${AMICI_DIR}/include
    ${AMICI_DIR}
    ../loadbalancer/src/
    ../common/src
    ../optimization/src
    ${SUITESPARSE_INCLUDE_DIRS}
    ${SUNDIALS_INCLUDE_DIRS}

)

add_subdirectory(examples)
add_subdirectory(src)

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif(${BUILD_TESTS})
