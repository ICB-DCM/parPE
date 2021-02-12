project(parpeoptimization)

find_package(Ceres COMPONENTS
    HINTS "${CMAKE_SOURCE_DIR}/ThirdParty/ceres-solver-2.0.0/build/install")

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

add_library(${PROJECT_NAME} ${SRC_LIST})

if(${PARPE_ENABLE_IPOPT})
    target_sources(${PROJECT_NAME}
        PRIVATE localOptimizationIpopt.cpp
        PRIVATE localOptimizationIpoptTNLP.cpp)
    list(APPEND HEADER_LIST
        localOptimizationIpopt.h
        localOptimizationIpoptTNLP.h)
    target_link_libraries(${PROJECT_NAME} PUBLIC PkgConfig::IPOPT)
endif(${PARPE_ENABLE_IPOPT})

if(${PARPE_ENABLE_CERES})
    target_sources(${PROJECT_NAME}
        PRIVATE localOptimizationCeres.cpp)
    list(APPEND HEADER_LIST localOptimizationCeres.h)
    target_link_libraries(${PROJECT_NAME}
        PUBLIC ceres)
    set(PARPE_CERES_MINIGLOG_REDIRECT FALSE CACHE BOOL
        "Custom treatment of Ceres output? Requires miniglog headers in ceres directory.")
    if(${PARPE_CERES_MINIGLOG_REDIRECT})
        target_compile_definitions(${PROJECT_NAME}
            PUBLIC PARPE_CERES_MINIGLOG_REDIRECT)
    endif(${PARPE_CERES_MINIGLOG_REDIRECT})
endif(${PARPE_ENABLE_CERES})

if(${PARPE_ENABLE_DLIB})
    target_sources(${PROJECT_NAME}
        PRIVATE localOptimizationDlib.cpp)
    list(APPEND HEADER_LIST localOptimizationDlib.h)
    find_package(dlib REQUIRED
        HINTS "${CMAKE_SOURCE_DIR}/ThirdParty/dlib-19.7/build/install/lib/cmake/dlib/")
    target_include_directories(${PROJECT_NAME}
        PUBLIC ${dlib_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME}
        PUBLIC ${dlib_LIBRARIES})
endif(${PARPE_ENABLE_DLIB})

if(${PARPE_ENABLE_TOMS611})
    target_sources(${PROJECT_NAME}
        PRIVATE localOptimizationToms611.cpp)
    list(APPEND HEADER_LIST localOptimizationToms611.h)
    add_subdirectory(${CMAKE_SOURCE_DIR}/ThirdParty/toms611 toms611)
    target_include_directories(${PROJECT_NAME}
        PUBLIC ${CMAKE_SOURCE_DIR}/ThirdParty/toms611/)
    target_link_libraries(${PROJECT_NAME}
        PUBLIC toms611)
endif(${PARPE_ENABLE_TOMS611})

if(${PARPE_ENABLE_FSQP})
    target_sources(${PROJECT_NAME}
        PRIVATE localOptimizationFsqp.cpp)
    list(APPEND HEADER_LIST localOptimizationFsqp.h)
    add_subdirectory(${CMAKE_SOURCE_DIR}/ThirdParty/ffsqp ffsqp)
    target_include_directories(${PROJECT_NAME}
        PUBLIC ${CMAKE_SOURCE_DIR}/ThirdParty/ffsqp/)
    target_link_libraries(${PROJECT_NAME}
        PUBLIC fsqp)
endif(${PARPE_ENABLE_FSQP})

target_include_directories(${PROJECT_NAME}
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${HDF5_INCLUDE_DIRS}
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC parpecommon
    PUBLIC ${HDF5_HL_LIBRARIES}
    PUBLIC ${HDF5_C_LIBRARIES}
    PUBLIC ${CMAKE_THREAD_LIBS_INIT}
)

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets ARCHIVE DESTINATION lib)

string(REGEX REPLACE "([^;]+)"
    "${CMAKE_SOURCE_DIR}/include/parpeoptimization/\\1"
    HEADER_LIST "${HEADER_LIST}")

install(
    FILES ${HEADER_LIST}
    DESTINATION
    include/parpeoptimization)
