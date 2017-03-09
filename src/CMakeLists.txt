project(pe)

set(SRC_LIST
    dataprovider.c
    localoptimizationipopt.c
    objectivefunction.c
    misc.c
    resultwriter.c
    queuemaster.c
    queueworker.c
    simulationworker.c
    masterthread.c
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
