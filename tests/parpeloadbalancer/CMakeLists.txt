project(unittests_loadbalancer)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

set(SRC_LIST
    loadBalancerMasterTest.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    parpeloadbalancer
    gmock_main
    ${GCOV_LIBRARY}
)
