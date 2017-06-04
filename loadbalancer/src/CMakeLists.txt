project(loadbalancer)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SRC_LIST
    loadBalancerMaster.c
    loadBalancerWorker.c
    queue.c
)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    Threads::Threads
)
