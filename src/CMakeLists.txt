project(pe)

set(SRC_LIST
    dataprovider.c
    localoptimizationipopt.c
    objectivefunction.c
    misc.c
    mpiworker.c
    resultwriter.c
    masterqueue.c
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
