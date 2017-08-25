project(parpeloadbalancer)

set(SRC_LIST
    LoadBalancerWorker.cpp
    LoadBalancerMaster.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${CMAKE_THREAD_LIBS_INIT}
)
