include_directories(src/)

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif(${BUILD_TESTS})

add_subdirectory(src/)
add_subdirectory(examples/loadbalancer)
