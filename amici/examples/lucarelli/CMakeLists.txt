# Build Lucarelli example

# This requires python to generate the model code from SBML
find_package(PythonInterp 3.6)

# Generate model code and data; build model
include(ExternalProject)
set(MODEL_NAME lucarelli_12)
ExternalProject_Add(${MODEL_NAME}
    DOWNLOAD_DIR        "${CMAKE_CURRENT_BINARY_DIR}/${MODEL_NAME}-prefix/src/${MODEL_NAME}"
    DOWNLOAD_COMMAND  ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/build_lucarelli_model.py"
    # Don't install anything
    INSTALL_COMMAND   ""
    TEST_COMMAND      ""
    #BUILD_ALWAYS      1
    DEPENDS parpeamici
    )

# Rebuild if amici files are updated
ExternalProject_Add_StepDependencies(${MODEL_NAME} build parpeamici)

ExternalProject_Get_property(${MODEL_NAME} SOURCE_DIR)
ExternalProject_Get_property(${MODEL_NAME} BINARY_DIR)
message("Generating Lucarelli example model and data in ${SOURCE_DIR}")

# Import external library
set(MODEL_LIBRARIES "model_${MODEL_NAME}")
add_library(${MODEL_LIBRARIES} STATIC IMPORTED)
set_target_properties(model_${MODEL_NAME}
    PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES
    "${SOURCE_DIR}"
    IMPORTED_LOCATION
    "${BINARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}${MODEL_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX}"
    )


# estimateLucarelli executable
project(estimateLucarelli)

set(SRC_LIST
    main.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST})
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
# /estimateLucarelli executable


# simulateLucarelli executable
project(simulateLucarelli)

set(SRC_LIST
    main_simulate.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST})
add_dependencies(${PROJECT_NAME} ${MODEL_NAME})

target_link_libraries(${PROJECT_NAME}
    ${MODEL_LIBRARIES}
    parpeamici
    ${GCOV_LIBRARY}
)
# /simulateLucarelli executable_simulator
