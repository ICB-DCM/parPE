project(parpeoptimization)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-write-strings")

set(SRC_LIST
    localOptimizationIpopt.cpp
    localOptimizationIpoptTNLP.cpp
    localOptimizationCeres.cpp
    optimizationProblem.cpp
    multiStartOptimization.cpp
    optimizationResultWriter.cpp
    optimizer.cpp
    optimizationOptions.cpp
)

if(${ENABLE_DLIB})
    set(SRC_LIST ${SRC_LIST} localOptimizationDlib.cpp)
endif(${ENABLE_DLIB})


add_library(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    parpecommon
    ${IPOPT_LIBRARIES}
    ${CERES_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)

if(${ENABLE_DLIB})
    set(DLIB_NO_GUI_SUPPORT ON)
    set(DLIB_GIF_SUPPORT OFF)
    set(DLIB_JPEG_SUPPORT OFF)
    set(DLIB_PNG_SUPPORT OFF)
    set(DLIB_LINK_WITH_SQLITE3 OFF)
    add_subdirectory(${DLIB_DIR} dlib_build)
    target_link_libraries(${PROJECT_NAME} dlib::dlib)
endif(${ENABLE_DLIB})

if(${ENABLE_TOMS611})
    target_link_libraries(${PROJECT_NAME} toms611)
endif(${ENABLE_TOMS611})


get_target_property(PARPEOPTIMIZATION_INCLUDE_DIRS ${PROJECT_NAME} INCLUDE_DIRECTORIES)
set_property(GLOBAL PROPERTY PARPEOPTIMIZATION_INCLUDE_DIRS "${PARPELOADBALANCER_INCLUDE_DIRS}")

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets ARCHIVE DESTINATION lib)
