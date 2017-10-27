include_directories(
    src/
    ../loadbalancer/src
    tests/
    ${IPOPT_INCLUDE_DIRS}
    ${CERES_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

if(${ENABLE_TOMS611})
    add_subdirectory(../ThirdParty/toms611/ toms611)
    include_directories(../ThirdParty/toms611/)

endif(${ENABLE_TOMS611})

add_subdirectory(src)

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif(${BUILD_TESTS})
