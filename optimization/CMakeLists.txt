include_directories(
    src/
    tests/
    ${IPOPT_INCLUDE_DIRS}
    ${CERES_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

add_subdirectory(src)

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif(${BUILD_TESTS})
