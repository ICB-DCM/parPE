project(parpeamici)

set(SRC_LIST
    MultiConditionDataProvider.cpp
    multiConditionProblem.cpp
    multiConditionProblemResultWriter.cpp
    simulationWorkerAmici.cpp
    steadystateSimulator.cpp
    amiciUserDataSerialization.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}

"${SUNDIALS_LIB_DIR}/libsundials_nvecserial${SHARED_OBJECT_EXTENSION}"
"${SUNDIALS_LIB_DIR}/libsundials_cvodes${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libcolamd${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libklu${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libbtf${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libamd${SHARED_OBJECT_EXTENSION}"
"${SUITESPARSE_LIB_DIR}/libsuitesparseconfig${SHARED_OBJECT_EXTENSION}"

    parpecommon
    ${IPOPT_LIBRARIES}
    ${CERES_LIBRARIES}
    ${SUNDIALS_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${Boost_SERIALIZATION_LIBRARY_RELEASE}
    ${CMAKE_THREAD_LIBS_INIT}
)
