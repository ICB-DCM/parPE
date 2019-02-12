find_package(Amici HINTS ${CMAKE_SOURCE_DIR}/deps/AMICI/build)

# create python virtual environment
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/venv
    COMMENT "Setting up virtual environment..."
    COMMAND sh -c "${CMAKE_SOURCE_DIR}/misc/venv.sh ${CMAKE_BINARY_DIR}")
# for manual invocation
add_custom_target(setup_venv
    COMMENT "Setting up virtual environment..."
    COMMAND sh -c "${CMAKE_SOURCE_DIR}/misc/venv.sh ${CMAKE_BINARY_DIR}")

add_subdirectory(steadystate)
add_subdirectory(lucarelli)
add_subdirectory(raia)
