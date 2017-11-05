project(example_loadbalancer)


set(SRC_LIST_EXE
    main.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE})

target_link_libraries(${PROJECT_NAME}
    ${CMAKE_THREAD_LIBS_INIT}
    parpeloadbalancer
)


if(${BUILD_TESTS})
    # --allow-run-as-root is required for running on shippable systems
    add_test (NAME ${PROJECT_NAME}
        COMMAND ${TESTS_MPIEXEC_COMMAND} ${PROJECT_NAME}
        )
endif()
