# Steady-state model example

# This requires python to generate the model code from SBML
find_package(PythonInterp 3.6)

# Generate model code and data; build model
include(ExternalProject)
set(MODEL_NAME steadystate_scaled)
ExternalProject_Add(
    ${MODEL_NAME}
    DOWNLOAD_DIR    "${CMAKE_CURRENT_BINARY_DIR}/${MODEL_NAME}-prefix/src/${MODEL_NAME}"
    DOWNLOAD_COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/create_steadystate_amici_model.py"
    SOURCE_SUBDIR    model_steadystate_scaled
    # Don't install anything
    INSTALL_COMMAND  ""
    TEST_COMMAND     ""
    #BUILD_ALWAYS    1
    DEPENDS          parpeamici
)

# Rebuild if amici files are updated
ExternalProject_Add_StepDependencies(${MODEL_NAME} download createSteadystateExampleSBML.py)
ExternalProject_Add_StepDependencies(${MODEL_NAME} download create_steadystate_amici_model.py)
ExternalProject_Add_StepDependencies(${MODEL_NAME} build parpeamici)

ExternalProject_Get_property(${MODEL_NAME} SOURCE_DIR)
ExternalProject_Get_property(${MODEL_NAME} BINARY_DIR)
message("Generating ${MODEL_NAME} example model and data in ${SOURCE_DIR} (${BINARY_DIR})")

# Import external library
set(MODEL_LIBRARIES "model_${MODEL_NAME}")
add_library(${MODEL_LIBRARIES} STATIC IMPORTED)
set_target_properties(model_${MODEL_NAME}
    PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES
    "${SOURCE_DIR}"
    IMPORTED_LOCATION
    "${BINARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}${MODEL_LIBRARIES}${CMAKE_STATIC_LIBRARY_SUFFIX}"
)
message("Imported lib ${MODEL_LIBRARIES}")


# example_steadystate executable
project(example_steadystate)

set(SRC_LIST_EXE
    main.cpp
    steadystateProblem.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE})
add_dependencies(${PROJECT_NAME} ${MODEL_NAME})
target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    Upstream::amici
    parpeoptimization
    parpeloadbalancer
    ${CMAKE_THREAD_LIBS_INIT}
)
# /example_steadystate executable

# example_steadystate_parallel executable
if(${PARPE_ENABLE_MPI})
    project(example_steadystate_parallel)

    set(SRC_LIST_EXE_PARALLEL
        main_parallel.cpp
        steadystateProblem.cpp
        steadystateProblemParallel.cpp
        steadyStateMultiConditionDataprovider.cpp
        )

    add_executable(${PROJECT_NAME} ${SRC_LIST_EXE_PARALLEL})
    add_dependencies(${PROJECT_NAME} ${MODEL_NAME})
    target_link_libraries(${PROJECT_NAME}
        ${MODEL_LIBRARIES}
        parpeoptimization
        parpeloadbalancer
        parpeamici
        )
endif()
# /example_steadystate_parallel executable


# example_steadystate_multi executable
project(example_steadystate_multi)

set(SRC_LIST_EXE_PARALLEL
    main_multicondition.cpp
    steadyStateMultiConditionDataprovider.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE_PARALLEL})
add_dependencies(${PROJECT_NAME} ${MODEL_NAME})
target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    parpeamici
    ${GCOV_LIBRARY}
)

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
    )

if(${BUILD_TESTS} AND ${PARPE_ENABLE_MPI})
    add_test (NAME ${PROJECT_NAME}
        COMMAND ${TESTS_MPIEXEC_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME} -o ${PROJECT_NAME}-testresults ${SOURCE_DIR}/example_data.h5
        )
    add_test (NAME ${PROJECT_NAME}_gradient_check
        COMMAND ${TESTS_MPIEXEC_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME} -t gradient_check -o ${PROJECT_NAME}-testresults ${SOURCE_DIR}/example_data.h5
        )
endif()
# /example_steadystate_multi executable


# example_steadystate_multi_simulator executable

project(example_steadystate_multi_simulator)

set(SRC_LIST_EXE_PARALLEL
    main_simulator.cpp
    steadyStateMultiConditionDataprovider.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE_PARALLEL})
add_dependencies(${PROJECT_NAME} ${MODEL_NAME})
target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    parpeamici
    ${GCOV_LIBRARY}
)


if(${BUILD_TESTS})
#    add_test (NAME ${PROJECT_NAME}
#        COMMAND ${TESTS_MPIEXEC_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME} -o ${PROJECT_NAME}-testresults ${CMAKE_CURRENT_LIST_DIR}/data.h5
#        )
endif()
# /example_steadystate_multi executable_simulator


if(${BUILD_TESTS})
    project(test_steadystate)

    add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

    set(SRC_LIST_CPP
        ../../tests/main.cpp
        ../../../common/src/testingMisc.cpp
        steadystateProblem.cpp
        exampleSteadystateScaledTest.cpp
    )

    add_executable(${PROJECT_NAME} ${SRC_LIST_CPP})
    add_dependencies(${PROJECT_NAME} ${MODEL_NAME})
    target_link_libraries(${PROJECT_NAME}
        parpeamici
        ${MODEL_LIBRARIES}
        ${CppUTest_LIBRARIES}
        ${GCOV_LIBRARY}
    )
endif()
