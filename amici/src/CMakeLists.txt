project(parpeamici)

set(SRC_LIST
    multiConditionDataProvider.cpp
    multiConditionProblem.cpp
    multiConditionProblemResultWriter.cpp
    steadystateSimulator.cpp
    optimizationApplication.cpp
    simulationRunner.cpp
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
