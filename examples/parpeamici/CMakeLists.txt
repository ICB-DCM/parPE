find_package(Amici HINTS ${CMAKE_SOURCE_DIR}/deps/AMICI/build)

# create python virtual environment
add_custom_target(setup_venv)
add_custom_command(
    TARGET setup_venv PRE_BUILD
    COMMENT "Setting up virtual environment..."
    COMMAND sh -c "${CMAKE_SOURCE_DIR}/misc/venv.sh ${CMAKE_BINARY_DIR}")
# for manual invocation


add_subdirectory(steadystate)
add_subdirectory(lucarelli)
add_subdirectory(raia)
