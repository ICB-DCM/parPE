project(pe)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SRC_LIST
    localOptimizationIpopt.cpp
    localOptimizationCeres.cpp
    misc.c
    loadBalancerMaster.c
    loadBalancerWorker.c
    queue.c
    optimizationProblem.cpp
    multiStartOptimization.cpp
    logging.c
    hdf5Misc.c
)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    ${IPOPT_LIBRARIES}
    ${CERES_LIBRARIES}
    ${SUNDIALS_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    Threads::Threads
)

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif(${BUILD_TESTS})

add_subdirectory(examples/steadystate)
add_subdirectory(examples/loadbalancer)
