include_directories(${CMAKE_CURRENT_SOURCE_DIR})

project(unittests_optimization)

set(SRC_LIST_CPP
    ../parpecommon/testingMisc.cpp
    quadraticTestProblem.cpp
    multiStartOptimizationTest.cpp
    minibatchOptimizationTest.cpp
    optimizationResultWriterTest.cpp
    optimizationOptionsTest.cpp
    optimizationProblemTest.cpp)

add_executable(${PROJECT_NAME} ${SRC_LIST_CPP})

target_link_libraries(${PROJECT_NAME} parpeoptimization parpecommon gmock_main
                      ${GCOV_LIBRARY})

if(${PARPE_ENABLE_FIDES})
  target_sources(${PROJECT_NAME} PRIVATE localOptimizationFidesTest.cpp)
endif()

if(${PARPE_ENABLE_IPOPT})
  target_sources(${PROJECT_NAME} PRIVATE localOptimizationIpoptTest.cpp)
  target_link_libraries(${PROJECT_NAME} PkgConfig::IPOPT)
endif()

if(${PARPE_ENABLE_CERES})
  target_sources(${PROJECT_NAME} PRIVATE localOptimizationCeresTest.cpp)
endif()

if(${PARPE_ENABLE_TOMS611})
  target_sources(${PROJECT_NAME} PRIVATE localOptimizationToms611Test.cpp)
  target_link_libraries(${PROJECT_NAME} toms611)
endif(${PARPE_ENABLE_TOMS611})

if(${PARPE_ENABLE_FSQP})
  target_sources(${PROJECT_NAME} PRIVATE localOptimizationFsqpTest.cpp)
  target_link_libraries(${PROJECT_NAME} fsqp)
endif(${PARPE_ENABLE_FSQP})

gtest_discover_tests(${PROJECT_NAME})
