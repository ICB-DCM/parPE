project(example_loadbalancer)

set(SRC_LIST_EXE
    main.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE})

target_link_libraries(${PROJECT_NAME}
    Threads::Threads
    pe
)
