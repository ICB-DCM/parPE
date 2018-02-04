project(parpeamici)

set(SRC_LIST
    MultiConditionDataProvider.cpp
    multiConditionProblem.cpp
    multiConditionProblemResultWriter.cpp
    simulationWorkerAmici.cpp
    steadystateSimulator.cpp
    optimizationApplication.cpp
    SimulationRunner.cpp
    simulationResultWriter.cpp
    standaloneSimulator.cpp
    amiciMisc.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    parpeoptimization
    parpeloadbalancer
    parpecommon
    Upstream::amici
    ${Boost_SERIALIZATION_LIBRARY_RELEASE}
    ${CMAKE_THREAD_LIBS_INIT}
)

get_target_property(PARPEAMICI_INCLUDE_DIRS ${PROJECT_NAME} INCLUDE_DIRECTORIES)
set_property(GLOBAL PROPERTY PARPEAMICI_INCLUDE_DIRS "${PARPEAMICI_INCLUDE_DIRS}")

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)
