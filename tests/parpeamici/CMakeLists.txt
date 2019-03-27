project(unittests_amici)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

set(SRC_LIST_CPP
    main.cpp
    ../parpecommon/testingMisc.cpp
    ../parpeoptimization/quadraticTestProblem.cpp
    amiciSimulationRunnerTest.h
    multiConditionDataProviderTest.h
    multiConditionProblemTest.h
    simulationResultWriterTest.h
    hierarchicalOptimizationTest.h
    ${GTestSrc}/src/gtest-all.cc
    ${GMockSrc}/src/gmock-all.cc
)

add_executable(${PROJECT_NAME} ${SRC_LIST_CPP})

# generate test h5 file first
add_custom_target(prepare_test_hierarchical_optimization
    COMMAND ${CMAKE_SOURCE_DIR}/misc/run_in_venv.sh ${CMAKE_BINARY_DIR}/venv ${CMAKE_CURRENT_SOURCE_DIR}/hierarchicalOptimizationTest.py
    COMMENT "Creating test data using hierarchicalOptimizationTest.py"
    DEPENDS setup_venv)
add_dependencies(${PROJECT_NAME} prepare_test_hierarchical_optimization)

#target_compile_options(${PROJECT_NAME} PUBLIC
#    -include ${CMAKE_CURRENT_LIST_DIR}/../../common/src/STLCompatibleMemoryLeakDetectorMacros.h)

target_link_libraries(${PROJECT_NAME}
    parpeamici
    ${GCOV_LIBRARY}
)
