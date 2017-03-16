project(benchmarkProblem)

set(SRC_LIST
    ${AMICI_DIR}/src/ami_hdf5.cpp
    dataprovider.c
    objectiveFunction.c
    problem.c
    resultwriter.c
    simulationWorker.c
    masterthread.c
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
    pe
)



project(cpp_benchmark)

set(SRC_LIST_EXE
    main.c
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE})

target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    ${IPOPT_LIBRARIES}
    ${SUNDIALS_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${MPE_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    -lpthread
    benchmarkProblem
    benchmark_model_simulation_lib
    pe
)
