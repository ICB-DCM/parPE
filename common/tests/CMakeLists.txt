project(unittests_common)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${CMAKE_CURRENT_LIST_DIR}/../../common/src/STLCompatibleMemoryLeakDetectorMacros.h")

include_directories(${CPPUTEST_INCLUDE_DIR})

set(SRC_LIST
    main.cpp
    commonTests.cpp
    hdf5MiscTests.cpp
    ../src/testingMisc.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${CPPUTEST_LIBRARY}
    ${CMAKE_THREAD_LIBS_INIT}
    parpecommon
    ${GCOV_LIBRARY}
)
