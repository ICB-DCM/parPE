project(unittests_pe)

add_test (NAME ${PROJECT_NAME}_c COMMAND ${PROJECT_NAME}_c -c)
add_test (NAME ${PROJECT_NAME}_cpp COMMAND ${PROJECT_NAME}_cpp -c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -include ${CPPUTEST_DIR}/include/CppUTest/MemoryLeakDetectorMallocMacros.h")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${CPPUTEST_DIR}/include/CppUTest/MemoryLeakDetectorNewMacros.h")

include_directories(${CPPUTEST_INCLUDE_DIR})




set(SRC_LIST_CPP
    main.cpp
    ../testingMisc.c
    ../localOptimizationIpoptTest.cpp
    ../localOptimizationCeresTest.cpp
    ../multiStartOptimizationTest.cpp
    quadraticTestProblem.cpp
)

add_executable(${PROJECT_NAME}_cpp ${SRC_LIST_CPP})

target_link_libraries(${PROJECT_NAME}_cpp
    pe
    ${IPOPT_LIBRARIES}
    ${CERES_LIBRARIES}
    ${SUNDIALS_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${CPPUTEST_LIBRARY}
    -lpthread
    ${BLAS_LIBRARIES}
)
