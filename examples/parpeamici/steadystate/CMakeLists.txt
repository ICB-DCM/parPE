# Steady-state model example

# This requires python to generate the model code from SBML
find_package(Python3 COMPONENTS Interpreter)
# TODO: check python version

# Generate model code and data; build model
include(ExternalProject)
set(MODEL_NAME steadystate_scaled)
set(BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${MODEL_NAME}-prefix/src/${MODEL_NAME}-build)
set(MODEL_LIBRARIES "model_${MODEL_NAME}")
set(MODEL_LIBRARY "${BINARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}${MODEL_LIBRARIES}${CMAKE_STATIC_LIBRARY_SUFFIX}")
ExternalProject_Add(
    ${MODEL_NAME}
    DOWNLOAD_DIR    "${CMAKE_CURRENT_BINARY_DIR}/${MODEL_NAME}-prefix/src/${MODEL_NAME}"
    DOWNLOAD_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/recreateExampleData.sh ${CMAKE_BINARY_DIR}
    SOURCE_SUBDIR    model_steadystate_scaled
    CMAKE_ARGS    -DAmici_DIR=${Amici_DIR}
    # Don't install anything
    INSTALL_COMMAND  ""
    TEST_COMMAND     ""
    #BUILD_ALWAYS    1
    DEPENDS          parpeamici
    BUILD_BYPRODUCTS "${MODEL_LIBRARY}"
)

# Rebuild if amici files are updated
ExternalProject_Add_StepDependencies(${MODEL_NAME} download
    recreateExampleData.sh
    ${CMAKE_SOURCE_DIR}/python/parpe/hdf5_pe_input.py
    ${CMAKE_SOURCE_DIR}/examples/parpeamici/steadystate/create_steadystate_amici_model.py
    ${CMAKE_SOURCE_DIR}/examples/parpeamici/steadystate/createSteadystateExampleSBML.py)
ExternalProject_Add_StepDependencies(${MODEL_NAME} build parpeamici)
ExternalProject_Add_StepDependencies(${MODEL_NAME} download setup_venv)

ExternalProject_Get_property(${MODEL_NAME} SOURCE_DIR)
ExternalProject_Get_property(${MODEL_NAME} BINARY_DIR)
message(STATUS "Generating ${MODEL_NAME} example model and data in ${SOURCE_DIR} (${BINARY_DIR})")

# Import external library
add_library(${MODEL_LIBRARIES} STATIC IMPORTED)
# Work-around for 'Imported target includes non-existent path'
file(MAKE_DIRECTORY "${SOURCE_DIR}/model_steadystate_scaled")
set_target_properties(model_${MODEL_NAME}
    PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES
    "${SOURCE_DIR}/model_steadystate_scaled"
    IMPORTED_LOCATION "${MODEL_LIBRARY}"
)
message(STATUS "Imported lib ${MODEL_LIBRARIES}")


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
)
# /example_steadystate executable


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

if("${BUILD_TESTING}" AND ${PARPE_ENABLE_MPI})
    add_test (NAME pytest_steadystate_example
        COMMAND ${CMAKE_SOURCE_DIR}/misc/run_in_venv.sh
            ${CMAKE_BINARY_DIR}/venv
            pytest -vv ${CMAKE_CURRENT_LIST_DIR}/test_steadystate.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
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


if(${BUILD_TESTING})
#    add_test (NAME ${PROJECT_NAME}
#        COMMAND ${TESTS_MPIEXEC_COMMAND} ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME} -o ${PROJECT_NAME}-testresults ${CMAKE_CURRENT_LIST_DIR}/data.h5
#        )
endif()
# /example_steadystate_multi executable_simulator


if(${BUILD_TESTING})
    project(test_steadystate)

    add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

    set(SRC_LIST_CPP
        exampleSteadystateScaledTest.h
        ${CMAKE_SOURCE_DIR}/tests/parpecommon/testingMisc.cpp
        steadystateProblem.cpp
    )

    add_executable(${PROJECT_NAME} ${SRC_LIST_CPP})
    add_dependencies(${PROJECT_NAME} ${MODEL_NAME})
    target_link_libraries(${PROJECT_NAME}
        parpeamici
        gmock_main
        ${MODEL_LIBRARIES}
        ${GCOV_LIBRARY}
    )
endif()
