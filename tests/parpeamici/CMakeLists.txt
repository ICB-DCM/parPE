project(unittests_amici)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)



set(SRC_LIST_CPP
    main.cpp
    amiciSimulationRunnerTest.cpp
    multiConditionDataProviderTest.cpp
    multiConditionProblemTest.cpp
    simulationResultWriterTest.cpp
    hierarchicalOptimizationTest.cpp
    ../parpecommon/testingMisc.cpp
    ../parpeoptimization/quadraticTestProblem.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_CPP})

# generate test h5 file first
add_custom_target(prepare_test_hierarchical_optimization
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/hierarchicalOptimizationTest.py
    COMMENT "Creating test data using hierarchicalOptimizationTest.py")
add_dependencies(${PROJECT_NAME} prepare_test_hierarchical_optimization)

#target_compile_options(${PROJECT_NAME} PUBLIC
#    -include ${CMAKE_CURRENT_LIST_DIR}/../../common/src/STLCompatibleMemoryLeakDetectorMacros.h)

target_link_libraries(${PROJECT_NAME}
    parpeamici
    ${CppUTest_LIBRARIES}
    ${GCOV_LIBRARY}
)
