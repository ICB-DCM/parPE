project(unittests_amici)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

# Cannot use with boost serialization
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -include ${CPPUTEST_DIR}/include/CppUTest/MemoryLeakDetectorMallocMacros.h")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${CPPUTEST_DIR}/include/CppUTest/MemoryLeakDetectorNewMacros.h")

set(SRC_LIST_CPP
    main.cpp
    ../src/amiciUserDataSerializationTest.cpp
    ../src/simulationWorkerAmiciTest.cpp
    ../src/multiConditionProblemResultWriterTest.cpp
    ../src/multiConditionDataProviderTest.cpp
    ../src/multiConditionProblemTest.cpp
    ../../common/src/testingMisc.c
)
set(SHARED_OBJECT_EXTENSION .so)

add_executable(${PROJECT_NAME} ${SRC_LIST_CPP})

target_link_libraries(${PROJECT_NAME}
    parpeamici
    parpeoptimization
    parpeloadbalancer
    parpecommon
    ${AMICI_LIBRARIES}
    ${IPOPT_LIBRARIES}
    ${CERES_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${CPPUTEST_LIBRARY}
    ${BLAS_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${GCOV_LIBRARY}
)
