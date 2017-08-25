project(unittests_loadbalancer)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -include ${CPPUTEST_DIR}/include/CppUTest/MemoryLeakDetectorMallocMacros.h")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include /usr/include/c++/6/bits/stl_tree.h -include ${CPPUTEST_DIR}/include/CppUTest/MemoryLeakDetectorNewMacros.h")

include_directories(${CPPUTEST_INCLUDE_DIR})

set(SRC_LIST_C
    main_c.cpp
    ../src/loadBalancerWorkerTest.c
    ../src/queueTest.c
)

set(SRC_LIST_CPP
    main.cpp
#    ../src/loadBalancerWorkerTest.cpp
    ../src/loadBalancerMasterTest.cpp

)

add_executable(${PROJECT_NAME} ${SRC_LIST_CPP})

target_link_libraries(${PROJECT_NAME}
    ${CPPUTEST_LIBRARY}
    ${CMAKE_THREAD_LIBS_INIT}
    parpeloadbalancer
    ${GCOV_LIBRARY}
)
