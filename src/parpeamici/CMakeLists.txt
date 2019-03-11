find_package(Amici HINTS ${CMAKE_SOURCE_DIR}/deps/AMICI/build)

if(CMAKE_BUILD_TYPE MATCHES Release)
    set(Boost_USE_STATIC_LIBS TRUE)
else()
    # For python module we need -fPIC which is only the case with shared libs
    set(Boost_USE_STATIC_LIBS FALSE)
endif()
find_package(Boost COMPONENTS serialization REQUIRED)

project(parpeamici)

set(SRC_LIST
    multiConditionDataProvider.cpp
    multiConditionProblem.cpp
    steadystateSimulator.cpp
    optimizationApplication.cpp
    amiciSimulationRunner.cpp
    simulationResultWriter.cpp
    standaloneSimulator.cpp
    amiciMisc.cpp
    hierarchicalOptimization.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_include_directories(${PROJECT_NAME}
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${Boost_INCLUDE_DIR}
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC parpeoptimization
    PUBLIC parpeloadbalancer
    PUBLIC parpecommon
    PUBLIC Upstream::amici
    PUBLIC ${Boost_SERIALIZATION_LIBRARY_RELEASE}
)

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)

install(
    DIRECTORY
    ${CMAKE_SOURCE_DIR}/include/parpeamici
    DESTINATION
    include)
