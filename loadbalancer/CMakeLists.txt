add_subdirectory(src/)

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif(${BUILD_TESTS})

if(${BUILD_EXAMPLES})
    add_subdirectory(examples/loadbalancer)
endif(${BUILD_EXAMPLES})
