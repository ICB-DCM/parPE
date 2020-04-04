project(example_loadbalancer)


set(SRC_LIST_EXE
    main.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE})

target_link_libraries(${PROJECT_NAME}
    ${CMAKE_THREAD_LIBS_INIT}
    parpeloadbalancer
)


if(${BUILD_TESTING})
    add_test (NAME ${PROJECT_NAME}
        COMMAND ${TESTS_MPIEXEC_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}
        )
endif()
