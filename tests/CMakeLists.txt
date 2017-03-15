project(unittests)

add_test (NAME test1 COMMAND unittests)

set(CPPUTEST_DIR /home/dweindl/src/cpputest/)
set(CPPUTEST_INCLUDE_DIR ${CPPUTEST_DIR}/include)
set(CPPUTEST_LIBRARY ${CPPUTEST_DIR}/build/src/CppUTest/libCppUTest.a ${CPPUTEST_DIR}/build/src/CppUTestExt/libCppUTestExt.a)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -include ${CPPUTEST_DIR}/include/CppUTest/MemoryLeakDetectorMallocMacros.h")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${CPPUTEST_DIR}/include/CppUTest/MemoryLeakDetectorNewMacros.h")

set(SRC_LIST
    main.cpp
    misc.c
    objectiveFunctionTest.c
    testqueueworker.c
    testqueuemaster.c
)

include_directories(${CPPUTEST_INCLUDE_DIR})

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    ${IPOPT_LIBRARIES}
    ${SUNDIALS_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${MPE_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${CPPUTEST_LIBRARY}
    -lpthread
    pe
)
