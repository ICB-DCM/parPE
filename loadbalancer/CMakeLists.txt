add_subdirectory(src/)

if(${BUILD_TESTS} AND ${PARPE_ENABLE_MPI})
    add_subdirectory(tests)
endif()

if(${BUILD_EXAMPLES} AND ${PARPE_ENABLE_MPI})
    add_subdirectory(examples/loadbalancer)
endif()
