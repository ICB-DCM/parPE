project(unittests_optimization)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${CMAKE_CURRENT_LIST_DIR}/../../common/src/STLCompatibleMemoryLeakDetectorMacros.h")

set(SRC_LIST_CPP
    main.cpp
    ../src/localOptimizationIpoptTest.cpp
    ../src/localOptimizationCeresTest.cpp
    ../src/multiStartOptimizationTest.cpp
    ../src/optimizationResultWriterTest.cpp
    ../src/optimizationOptionsTest.cpp
    quadraticTestProblem.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_CPP})

target_link_libraries(${PROJECT_NAME}
    parpeoptimization
    parpecommon
    ${CPPUTEST_LIBRARY}
    ${GCOV_LIBRARY}
)

if(${ENABLE_TOMS611})
    project(unittests_optimization_toms611)

    add_executable(${PROJECT_NAME}
        main.cpp
        ../src/localOptimizationToms611Test.cpp
    )

    target_link_libraries(${PROJECT_NAME}
        toms611
        parpeoptimization
        parpecommon
        ${CPPUTEST_LIBRARY}
        ${GCOV_LIBRARY}
    )
endif(${ENABLE_TOMS611})
