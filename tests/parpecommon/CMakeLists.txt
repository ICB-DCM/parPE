project(unittests_common)

set(SRC_LIST
    main.cpp
    commonTests.h
    hdf5MiscTests.h
    ${CMAKE_SOURCE_DIR}/tests/parpecommon/testingMisc.cpp
    ${GTestSrc}/src/gtest-all.cc
    ${GMockSrc}/src/gmock-all.cc
)

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${CMAKE_THREAD_LIBS_INIT}
    parpecommon
    ${GCOV_LIBRARY}
)

gtest_discover_tests(${PROJECT_NAME})
