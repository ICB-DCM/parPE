project(unittests_loadbalancer)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

include_directories(${CPPUTEST_INCLUDE_DIRS})

set(SRC_LIST
    main.cpp
    loadBalancerMasterTest.h
    ${GTestSrc}/src/gtest-all.cc
    ${GMockSrc}/src/gmock-all.cc
)

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${CMAKE_THREAD_LIBS_INIT}
    parpeloadbalancer
    ${GCOV_LIBRARY}
)
