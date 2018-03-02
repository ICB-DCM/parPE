set(MODEL_LIBRARIES "model_steadystate")
set(MODEL_INCLUDE_DIRS "model")

include_directories(${MODEL_INCLUDE_DIRS})

# AMICI-generated model files
add_subdirectory(model)

# example_steadystate executable
project(example_steadystate)

set(SRC_LIST_EXE
    main.cpp
    steadystateProblem.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE})

target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    Upstream::amici
    parpeoptimization
    parpeloadbalancer
    ${CMAKE_THREAD_LIBS_INIT}
)
# /example_steadystate executable

# example_steadystate_parallel executable
project(example_steadystate_parallel)

set(SRC_LIST_EXE_PARALLEL
    main_parallel.cpp
    steadystateProblem.cpp
    steadystateProblemParallel.cpp
    steadyStateMultiConditionDataprovider.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE_PARALLEL})

target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    parpeoptimization
    parpeloadbalancer
    parpeamici
)
# /example_steadystate_parallel executable


# example_steadystate_multi executable

project(example_steadystate_multi)

set(SRC_LIST_EXE_PARALLEL
    main_multicondition.cpp
    steadyStateMultiConditionDataprovider.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE_PARALLEL})

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

if(${BUILD_TESTS})
    add_test (NAME ${PROJECT_NAME}
        COMMAND ${TESTS_MPIEXEC_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME} -o ${PROJECT_NAME}-testresults ${CMAKE_CURRENT_LIST_DIR}/data.h5
        )
    add_test (NAME ${PROJECT_NAME}_gradient_check
        COMMAND ${TESTS_MPIEXEC_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME} -t gradient_check -o ${PROJECT_NAME}-testresults ${CMAKE_CURRENT_LIST_DIR}/data.h5
        )
endif()
# /example_steadystate_multi executable


# example_steadystate_multi_simulator executable

project(example_steadystate_multi_simulator)

set(SRC_LIST_EXE_PARALLEL
    main_simulator.cpp
    steadyStateMultiConditionDataprovider.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE_PARALLEL})

target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    parpeamici
    ${GCOV_LIBRARY}
)


if(${BUILD_TESTS})
#    add_test (NAME ${PROJECT_NAME}
#        COMMAND ${TESTS_MPIEXEC_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME} -o ${PROJECT_NAME}-testresults ${CMAKE_CURRENT_LIST_DIR}/data.h5
#        )
endif()
# /example_steadystate_multi executable_simulator
