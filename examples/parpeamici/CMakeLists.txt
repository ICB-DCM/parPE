find_package(Amici HINTS ${CMAKE_SOURCE_DIR}/deps/AMICI/build)

add_subdirectory(steadystate)
add_subdirectory(lucarelli)
add_subdirectory(raia)
