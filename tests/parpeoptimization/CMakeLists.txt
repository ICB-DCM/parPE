include_directories(${CMAKE_CURRENT_SOURCE_DIR})

project(unittests_optimization)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

set(SRC_LIST_CPP
    main.cpp
    ../parpecommon/testingMisc.cpp
    quadraticTestProblem.cpp
    multiStartOptimizationTest.h
    minibatchOptimizationTest.h
    optimizationResultWriterTest.h
    optimizationOptionsTest.h
    optimizationProblemTest.h
    localOptimizationIpoptTest.h
    localOptimizationCeresTest.h
    ${GTestSrc}/src/gtest-all.cc
    ${GMockSrc}/src/gmock-all.cc
)

add_executable(${PROJECT_NAME} ${SRC_LIST_CPP})

target_link_libraries(${PROJECT_NAME}
    parpeoptimization
    parpecommon
    ${GCOV_LIBRARY}
)

if(${PARPE_ENABLE_TOMS611})
    target_link_libraries(${PROJECT_NAME} toms611)
endif(${PARPE_ENABLE_TOMS611})

if(${PARPE_ENABLE_FSQP})
    target_link_libraries(${PROJECT_NAME} fsqp)
endif(${PARPE_ENABLE_FSQP})
