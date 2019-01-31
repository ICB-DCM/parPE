project(unittests_loadbalancer)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

include_directories(${CPPUTEST_INCLUDE_DIRS})

set(SRC_LIST
    main.cpp
    loadBalancerMasterTest.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST})

#target_compile_options(${PROJECT_NAME} PUBLIC
#    -include ${CMAKE_CURRENT_LIST_DIR}/../../common/src/STLCompatibleMemoryLeakDetectorMacros.h)

target_link_libraries(${PROJECT_NAME}
    ${CppUTest_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    parpeloadbalancer
    ${GCOV_LIBRARY}
)
