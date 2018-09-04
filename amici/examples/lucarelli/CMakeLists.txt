set(MODEL_LIBRARIES "lucarelli_12")

# AMICI-generated model files
add_subdirectory(lucarelli_12)

# estimateLucarelli executable
project(estimateLucarelli)

set(SRC_LIST
    main.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    parpeamici
    ${GCOV_LIBRARY}
)

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
    )
# /estimateLucarelli executable


# simulateLucarelli executable
project(simulateLucarelli)

set(SRC_LIST
    main_simulate.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    parpeamici
    ${GCOV_LIBRARY}
)
# /simulateLucarelli executable_simulator
