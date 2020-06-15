project(parpeoptimization)

find_package(Ceres COMPONENTS
    HINTS "${CMAKE_SOURCE_DIR}/ThirdParty/ceres-solver-1.13.0/build/install")

if(${PARPE_ENABLE_CERES} AND NOT ${CERES_FOUND})
    message(FATAL_ERROR "Supposed to build with CERES support, "
        "but CERES not found.")
endif()

set(SRC_LIST
    optimizationProblem.cpp
    multiStartOptimization.cpp
    optimizationResultWriter.cpp
    optimizationOptions.cpp
    minibatchOptimization.cpp
)

set(HEADER_LIST
    minibatchOptimization.h
    multiStartOptimization.h
    optimizationOptions.h
    optimizationProblem.h
    optimizationResultWriter.h
    optimizer.h
    )

if(${PARPE_ENABLE_IPOPT})
    list(APPEND SRC_LIST localOptimizationIpopt.cpp localOptimizationIpoptTNLP.cpp)
    list(APPEND HEADER_LIST localOptimizationIpopt.h localOptimizationIpoptTNLP.h)
endif(${PARPE_ENABLE_IPOPT})

if(${PARPE_ENABLE_CERES})
    list(APPEND SRC_LIST localOptimizationCeres.cpp)
    list(APPEND HEADER_LIST localOptimizationCeres.h)
endif(${PARPE_ENABLE_CERES})

if(${PARPE_ENABLE_DLIB})
    list(APPEND SRC_LIST localOptimizationDlib.cpp)
    list(APPEND HEADER_LIST localOptimizationDlib.h)
endif(${PARPE_ENABLE_DLIB})

if(${PARPE_ENABLE_TOMS611})
    list(APPEND SRC_LIST localOptimizationToms611.cpp)
    list(APPEND HEADER_LIST localOptimizationToms611.h)
endif(${PARPE_ENABLE_TOMS611})

if(${PARPE_ENABLE_FSQP})
    list(APPEND SRC_LIST localOptimizationFsqp.cpp)
    list(APPEND HEADER_LIST localOptimizationFsqp.h)
endif(${PARPE_ENABLE_FSQP})

add_library(${PROJECT_NAME} ${SRC_LIST})

target_include_directories(${PROJECT_NAME}
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${HDF5_INCLUDE_DIRS}
)

if(${PARPE_ENABLE_CERES})
    target_include_directories(${PROJECT_NAME} PUBLIC ${CERES_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${CERES_LIBRARIES})
endif(${PARPE_ENABLE_CERES})

if(${PARPE_ENABLE_IPOPT})
    target_include_directories(${PROJECT_NAME} PUBLIC ${IPOPT_INCLUDE_DIRS})
endif(${PARPE_ENABLE_IPOPT})

set(PARPE_CERES_MINIGLOG_REDIRECT FALSE CACHE BOOL "Custom treatment of Ceres output? Requires miniglog headers in ceres directory.")
if(${PARPE_CERES_MINIGLOG_REDIRECT})
    target_compile_definitions(${PROJECT_NAME} PUBLIC PARPE_CERES_MINIGLOG_REDIRECT)
endif(${PARPE_CERES_MINIGLOG_REDIRECT})

if(${PARPE_ENABLE_TOMS611})
    add_subdirectory(${CMAKE_SOURCE_DIR}/ThirdParty/toms611 toms611)
    target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/ThirdParty/toms611/)
endif(${PARPE_ENABLE_TOMS611})

if(${PARPE_ENABLE_DLIB})
    find_package(dlib REQUIRED
        HINTS "${CMAKE_SOURCE_DIR}/ThirdParty/dlib-19.7/build/install/lib/cmake/dlib/")
    target_include_directories(${PROJECT_NAME} PUBLIC ${dlib_INCLUDE_DIRS})
endif(${PARPE_ENABLE_DLIB})

if(${PARPE_ENABLE_FSQP})
    add_subdirectory(${CMAKE_SOURCE_DIR}/ThirdParty/ffsqp ffsqp)
    target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/ThirdParty/ffsqp/)
endif(${PARPE_ENABLE_FSQP})


target_link_libraries(${PROJECT_NAME}
    parpecommon
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)

if(${PARPE_ENABLE_IPOPT})
    target_link_libraries(${PROJECT_NAME} ${IPOPT_LIBRARIES})
endif(${PARPE_ENABLE_IPOPT})

if(${PARPE_ENABLE_DLIB})
    target_link_libraries(${PROJECT_NAME} ${dlib_LIBRARIES})
endif(${PARPE_ENABLE_DLIB})

if(${PARPE_ENABLE_TOMS611})
    target_link_libraries(${PROJECT_NAME} toms611)
endif(${PARPE_ENABLE_TOMS611})

if(${PARPE_ENABLE_FSQP})
    target_link_libraries(${PROJECT_NAME} fsqp)
endif(${PARPE_ENABLE_FSQP})

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets ARCHIVE DESTINATION lib)

string(REGEX REPLACE "([^;]+)" "${CMAKE_SOURCE_DIR}/include/parpeoptimization/\\1" HEADER_LIST "${HEADER_LIST}")
install(
    FILES ${HEADER_LIST}
    DESTINATION
    include)
