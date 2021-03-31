project(unittests_common)

set(SRC_LIST
    commonTests.cpp
    hdf5MiscTests.cpp
    ${CMAKE_SOURCE_DIR}/tests/parpecommon/testingMisc.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${CMAKE_THREAD_LIBS_INIT}
    parpecommon
    gtest_main
    ${GCOV_LIBRARY}
)

gtest_discover_tests(${PROJECT_NAME})
