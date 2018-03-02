project(parpeoptimization)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-write-strings")

set(SRC_LIST
    localOptimizationIpopt.cpp
    localOptimizationIpoptTNLP.cpp
    localOptimizationCeres.cpp
    optimizationProblem.cpp
    multiStartOptimization.cpp
    optimizationResultWriter.cpp
    optimizationOptions.cpp
    minibatchOptimization.cpp
)

if(${PARPE_ENABLE_DLIB})
    set(SRC_LIST ${SRC_LIST} localOptimizationDlib.cpp)
endif(${PARPE_ENABLE_DLIB})

if(${PARPE_ENABLE_TOMS611})
    set(SRC_LIST ${SRC_LIST} localOptimizationToms611.cpp)
endif(${PARPE_ENABLE_TOMS611})

if(${PARPE_ENABLE_FSQP})
    set(SRC_LIST ${SRC_LIST} localOptimizationFsqp.cpp)
endif(${PARPE_ENABLE_FSQP})

add_library(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    parpecommon
    ${IPOPT_LIBRARIES}
    ${CERES_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)

if(${PARPE_ENABLE_DLIB})
    target_link_libraries(${PROJECT_NAME} ${dlib_LIBRARIES})
endif(${PARPE_ENABLE_DLIB})

if(${PARPE_ENABLE_TOMS611})
    target_link_libraries(${PROJECT_NAME} toms611)
endif(${PARPE_ENABLE_TOMS611})

if(${PARPE_ENABLE_FSQP})
    target_link_libraries(${PROJECT_NAME} fsqp)
endif(${PARPE_ENABLE_FSQP})

get_target_property(PARPEOPTIMIZATION_INCLUDE_DIRS ${PROJECT_NAME} INCLUDE_DIRECTORIES)
set_property(GLOBAL PROPERTY PARPEOPTIMIZATION_INCLUDE_DIRS "${PARPELOADBALANCER_INCLUDE_DIRS}")

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets ARCHIVE DESTINATION lib)
