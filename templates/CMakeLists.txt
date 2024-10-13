# To show in IDE
add_custom_target(
  templates SOURCES CMakeLists.template.txt main.cpp main_debug.cpp
                    main_nominal.cpp main_simulate.cpp)

set_target_properties(
  templates
  PROPERTIES
    INCLUDE_DIRECTORIES
    "${CMAKE_SOURCE_DIR}/include/;${CMAKE_SOURCE_DIR}/deps/AMICI/include/")
