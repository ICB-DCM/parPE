project(unittests_optimization)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${CPPUTEST_DIR}/include/CppUTest/MemoryLeakDetectorMallocMacros.h -include ${CPPUTEST_DIR}/include/CppUTest/MemoryLeakDetectorNewMacros.h")

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
