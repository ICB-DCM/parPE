project(objectiveFunctionBenchmarkModel)

set(SRC_LIST
    ${AMICI_DIR}/src/ami_hdf5.cpp
    dataprovider.c
    objectiveFunction.c
)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    ${SUNDIALS_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${MPE_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    -lpthread
)
