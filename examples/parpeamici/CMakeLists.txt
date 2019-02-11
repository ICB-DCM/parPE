find_package(Amici HINTS ${CMAKE_SOURCE_DIR}/deps/AMICI/build)

# create python virtual environment
message(STATUS "Setting up virtual environment...")
execute_process(COMMAND sh -c "${CMAKE_SOURCE_DIR}/misc/venv.sh ${PROJECT_BINARY_DIR}")

add_subdirectory(steadystate)
add_subdirectory(lucarelli)
add_subdirectory(raia)
