find_package(SWIG REQUIRED 3)
find_package(PythonLibs 3.6 REQUIRED)
find_package(Amici HINTS ${CMAKE_SOURCE_DIR}/deps/AMICI/build)

include(${SWIG_USE_FILE})

include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${CMAKE_SOURCE_DIR}/include)

get_target_property(AMICI_INCLUDE_DIRS Upstream::amici
                    INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${AMICI_INCLUDE_DIRS})

set(SOURCES parpe.i)
set_property(SOURCE parpe.i PROPERTY CPLUSPLUS ON)

set(CMAKE_DEBUG_POSTFIX "")
swig_add_library(
  parpe
  TYPE MODULE
  LANGUAGE Python
  SOURCES ${SOURCES})

target_link_libraries(_parpe Upstream::amici parpeamici parpeoptimization
                      parpeloadbalancer parpecommon)
