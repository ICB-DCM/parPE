project(parpeoptimization)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-write-strings")

set(SRC_LIST
    localOptimizationIpopt.cpp
    localOptimizationIpoptTNLP.cpp
    localOptimizationCeres.cpp
    localOptimizationDlib.cpp
    optimizationProblem.cpp
    multiStartOptimization.cpp
    optimizationResultWriter.cpp
    optimizer.cpp
    optimizationOptions.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST})

add_subdirectory(${DLIB_DIR} dlib_build)


target_link_libraries(${PROJECT_NAME}
    parpecommon
    ${IPOPT_LIBRARIES}
    ${CERES_LIBRARIES}
    dlib::dlib
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)

get_target_property(PARPEOPTIMIZATION_INCLUDE_DIRS ${PROJECT_NAME} INCLUDE_DIRECTORIES)
set_property(GLOBAL PROPERTY PARPEOPTIMIZATION_INCLUDE_DIRS "${PARPELOADBALANCER_INCLUDE_DIRS}")

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)
