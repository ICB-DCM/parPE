project(parpeoptimization)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-write-strings")

set(SRC_LIST
    localOptimizationIpopt.cpp
    localOptimizationCeres.cpp
    optimizationProblem.cpp
    multiStartOptimization.cpp
    optimizationResultWriter.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    parpecommon
    ${IPOPT_LIBRARIES}
    ${CERES_LIBRARIES}
    ${SUNDIALS_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)
