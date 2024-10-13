# --- Add targets for cmake-format https://cmake-format.readthedocs.io/ ---

# Find all CMakeFiles files
set(ALL_CMAKE_FILES
    CMakeLists.txt
    examples/parpeloadbalancer/CMakeLists.txt
    examples/parpeamici/steadystate/CMakeLists.txt
    examples/parpeamici/CMakeLists.txt
    examples/CMakeLists.txt
    src/parpeloadbalancer/CMakeLists.txt
    src/parpeamici/CMakeLists.txt
    src/CMakeLists.txt
    src/parpecommon/CMakeLists.txt
    src/parpeoptimization/CMakeLists.txt
    swig/CMakeLists.txt
    templates/CMakeLists.txt
    tests/parpeloadbalancer/CMakeLists.txt
    tests/parpeamici/CMakeLists.txt
    tests/CMakeLists.txt
    tests/CMakeLists.txt.in
    tests/parpecommon/CMakeLists.txt
    tests/parpeoptimization/CMakeLists.txt
    cmake/clang-tools.cmake
    cmake/getVersion.cmake
    cmake/BuildOptimized.cmake
    cmake/BuildType.cmake
    cmake/split_version.cmake
    cmake/cmakelang-tools.cmake
    cmake/ConfigureVersion.cmake
    # cmake/CodeCoverage.cmake
    CTestConfig.cmake)

list(JOIN ALL_CMAKE_FILES " " ALL_CMAKE_FILES)

# --- cmake-format ---

# Try to find cmake-format and add target if successful
find_program(CMAKE_FORMAT "cmake-format")
set(CMAKE_FORMAT_OPTIONS "--tab-size 4")
if(CMAKE_FORMAT)
  add_custom_target(
    cmake-format
    COMMAND bash -c
            "${CMAKE_FORMAT} ${CMAKE_FORMAT_OPTIONS} -i ${ALL_CMAKE_FILES}"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running cmake-format")
else()
  message(STATUS "cmake-format was not found")
endif()

# --- cmake-lint ---

# Try to find cmake-lint and add target if successful
find_program(CMAKE_LINT "cmake-lint")
if(CMAKE_LINT)
  add_custom_target(
    cmake-lint
    COMMAND bash -c "${CMAKE_LINT} ${ALL_CMAKE_FILES}"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running cmake-lint")
else()
  message(STATUS "cmake-lint was not found")
endif()
