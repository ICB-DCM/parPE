project(unittests_amici)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

set(SRC_LIST_CPP
    main.cpp
    ../src/simulationRunnerTest.cpp
    ../src/multiConditionProblemResultWriterTest.cpp
    ../src/multiConditionDataProviderTest.cpp
    ../src/multiConditionProblemTest.cpp
    ../src/simulationResultWriterTest.cpp
    ../src/hierarchicalOptimizationTest.cpp
    ../../common/src/testingMisc.cpp
    ../../optimization/tests/quadraticTestProblem.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_CPP})

#target_compile_options(${PROJECT_NAME} PUBLIC
#    -include ${CMAKE_CURRENT_LIST_DIR}/../../common/src/STLCompatibleMemoryLeakDetectorMacros.h)

target_link_libraries(${PROJECT_NAME}
    parpeamici
    ${CppUTest_LIBRARIES}
    ${GCOV_LIBRARY}
)
