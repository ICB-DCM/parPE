find_package(Amici HINTS ${CMAKE_SOURCE_DIR}/deps/AMICI/build)

# create python virtual environment
message(STATUS "Setting up virtual environment...")
add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/venv
    COMMAND sh -c "${CMAKE_SOURCE_DIR}/misc/venv.sh ${PROJECT_BINARY_DIR}")
# for manual invocation
add_custom_target(setup_venv
        COMMAND sh -c "${CMAKE_SOURCE_DIR}/misc/venv.sh ${PROJECT_BINARY_DIR}")

add_subdirectory(steadystate)
add_subdirectory(lucarelli)
add_subdirectory(raia)
