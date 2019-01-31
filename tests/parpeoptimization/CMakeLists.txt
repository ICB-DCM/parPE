include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#General unit tests
project(unittests_optimization)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

set(SRC_LIST_CPP
    main.cpp
    multiStartOptimizationTest.cpp
    minibatchOptimizationTest.cpp
    optimizationResultWriterTest.cpp
    optimizationOptionsTest.cpp
    optimizationProblemTest.cpp
    ../parpecommon/testingMisc.cpp
    quadraticTestProblem.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_CPP})

#target_compile_options(${PROJECT_NAME} PUBLIC
#    -include ${CMAKE_CURRENT_LIST_DIR}/../../common/src/STLCompatibleMemoryLeakDetectorMacros.h)

target_link_libraries(${PROJECT_NAME}
    parpeoptimization
    parpecommon
    ${CppUTest_LIBRARIES}
    ${GCOV_LIBRARY}
)

#######################
# per-optimizer tests #
#######################

#######################
# IPOPT
#######################
if(${PARPE_ENABLE_IPOPT})
    project(unittests_optimization_ipopt)
    add_executable(${PROJECT_NAME}
        main.cpp
        localOptimizationIpoptTest.cpp
        quadraticTestProblem.cpp
        )
    target_link_libraries(${PROJECT_NAME}
        parpeoptimization
        parpecommon
        ${CppUTest_LIBRARIES}
        ${GCOV_LIBRARY}
        )
    add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)
endif(${PARPE_ENABLE_IPOPT})

#######################
# CERES
#######################
if(${PARPE_ENABLE_CERES})
    project(unittests_optimization_ceres)
    add_executable(${PROJECT_NAME}
        main.cpp
        localOptimizationCeresTest.cpp
        quadraticTestProblem.cpp
        )
    target_link_libraries(${PROJECT_NAME}
        parpeoptimization
        parpecommon
        ${CppUTest_LIBRARIES}
        ${GCOV_LIBRARY}
        )
    add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)
endif(${PARPE_ENABLE_CERES})

#######################
# TOMS611
#######################
if(${PARPE_ENABLE_TOMS611})
    project(unittests_optimization_toms611)
    add_executable(${PROJECT_NAME}
        main.cpp
        localOptimizationToms611Test.cpp
        quadraticTestProblem.cpp
    )
    target_link_libraries(${PROJECT_NAME}
        toms611
        parpeoptimization
        parpecommon
        ${CppUTest_LIBRARIES}
        ${GCOV_LIBRARY}
    )
    add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)
endif(${PARPE_ENABLE_TOMS611})

#######################
# FFSQP
#######################
if(${PARPE_ENABLE_FSQP})
    project(unittests_optimization_fsqp)
    add_executable(${PROJECT_NAME}
        main.cpp
        localOptimizationFsqpTest.cpp
        quadraticTestProblem.cpp
    )
    target_link_libraries(${PROJECT_NAME}
        fsqp
        parpeoptimization
        parpecommon
        ${CppUTest_LIBRARIES}
        ${GCOV_LIBRARY}
    )
    add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)
endif(${PARPE_ENABLE_FSQP})

#######################
