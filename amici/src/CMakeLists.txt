project(parpeamici)

set(SRC_LIST
    MultiConditionDataProvider.cpp
    multiConditionProblem.cpp
    multiConditionProblemResultWriter.cpp
    simulationWorkerAmici.cpp
    steadystateSimulator.cpp
    optimizationApplication.cpp
    #amiciUserDataSerialization.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    parpeoptimization
    parpeloadbalancer
    parpecommon
    amici
    ${Boost_SERIALIZATION_LIBRARY_RELEASE}
    ${CMAKE_THREAD_LIBS_INIT}
)
