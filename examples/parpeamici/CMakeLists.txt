find_package(Amici HINTS ${CMAKE_SOURCE_DIR}/deps/AMICI/build)

# rebuild amici
add_custom_target(rebuild_amici)
add_custom_command(
    TARGET rebuild_amici PRE_BUILD
    COMMENT "Rebuilding AMICI..."
    COMMAND sh -c "${CMAKE_SOURCE_DIR}/misc/rebuild_amici.sh ${CMAKE_BINARY_DIR}"
    SOURCES ${CMAKE_SOURCE_DIR}/deps/AMICI/.gitrepo)


# create python virtual environment
add_custom_target(setup_venv)
add_custom_command(
    TARGET setup_venv PRE_BUILD
    COMMENT "Setting up virtual environment..."
    COMMAND sh -c "${CMAKE_SOURCE_DIR}/misc/venv.sh ${CMAKE_BINARY_DIR}"
    SOURCES ${CMAKE_SOURCE_DIR}/misc/venv.sh)


if(${BUILD_TESTING})
    add_test (
        NAME "parpe_python_package_tests"
        COMMAND ${CMAKE_SOURCE_DIR}/misc/run_in_venv.sh
            ${CMAKE_BINARY_DIR}/venv
            pytest ${CMAKE_SOURCE_DIR}/python
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )
endif()


add_subdirectory(steadystate)
