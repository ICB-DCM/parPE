# Add targets for clang-format and clang-tidy ############

# Find all source files
execute_process(
  COMMAND
    sh -c
    "git ls-tree -r HEAD --name-only | grep -E '(\\.cpp$)|(\\.h$)' | tr '\n' ' '"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE ALL_CXX_SOURCE_FILES)

# ########### clang-format ############

# Try to find clang-format and add target if successful
find_program(CLANG_FORMAT "clang-format")
if(CLANG_FORMAT)
  add_custom_target(
    clang-format
    COMMAND bash -c "/usr/bin/clang-format -i ${ALL_CXX_SOURCE_FILES}"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
else()
  message(STATUS "clang-format was not found")
endif()

# ########### clang-tidy ############
# Try to find clang-tidy and add target if successful
find_program(CLANG_TIDY "clang-tidy")
if(CLANG_TIDY)
  add_custom_target(
    clang-tidy
    COMMAND
      sh -c
      "/usr/bin/clang-tidy amici/src/*.cpp -- -std=c++11 ${CMAKE_INCLUDE_PATH}"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
else()
  message(STATUS "clang-tidy was not found")
endif()
