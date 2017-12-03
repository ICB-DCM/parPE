project(unittests_amici)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

set(SRC_LIST_CPP
    main.cpp
    ../src/simulationWorkerAmiciTest.cpp
    ../src/multiConditionProblemResultWriterTest.cpp
    ../src/multiConditionDataProviderTest.cpp
    ../src/multiConditionProblemTest.cpp
    ../src/simulationResultWriterTest.cpp
    ../../common/src/testingMisc.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_CPP})

#target_compile_options(${PROJECT_NAME} PUBLIC
#    -include ${CMAKE_CURRENT_LIST_DIR}/../../common/src/STLCompatibleMemoryLeakDetectorMacros.h)

target_link_libraries(${PROJECT_NAME}
    parpeamici
    ${CPPUTEST_LIBRARY}
    ${GCOV_LIBRARY}
)
