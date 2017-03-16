project(pe)

set(SRC_LIST
    localOptimizationIpopt.c
    misc.c
    loadBalancerMaster.c
    loadBalancerWorker.c
    queue.c
    optimizationProblem.c
)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    ${IPOPT_LIBRARIES}
    ${SUNDIALS_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${MPE_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    -lpthread
)

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif(${BUILD_TESTS})
