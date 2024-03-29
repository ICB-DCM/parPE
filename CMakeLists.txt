cmake_minimum_required(VERSION 3.15) # pkg_search_module GLOBAL
cmake_policy(VERSION 3.15...3.28)

if(POLICY CMP0074)
  # Use package_ROOT environment variables
  cmake_policy(SET CMP0074 NEW)
endif(POLICY CMP0074)

project(parpe)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(BuildType) # Ensure CMAKE_BUILD_TYPE is always set
include(CTest)

set(CMAKE_DEBUG_POSTFIX "-dbg")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wno-unused-function")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function")

if(BUILD_PYTHON_MODULE)
    # Build PIC code to be used for swig/python module
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()


# <versioning>
# Version number from git
add_custom_target(get_version ALL
    COMMAND "${CMAKE_COMMAND}"
    "-DINFILE=${PROJECT_SOURCE_DIR}/src/parpecommon/parpeVersion.h.in"
    "-DOUTFILE=${CMAKE_BINARY_DIR}/src/parpecommon/parpeVersion.h"
    "-P" "${PROJECT_SOURCE_DIR}/cmake/ConfigureVersion.cmake")

include("getVersion")
include("split_version")
split_version(${PARPE_VERSION} libname major minor patch)
# </versioning>


# <Build options>
set(BUILD_EXAMPLES TRUE CACHE BOOL "Build examples?")
set(BUILD_DOCS FALSE CACHE BOOL "Generate documentation?")
set(BUILD_PYTHON_MODULE FALSE CACHE BOOL "Build Python module?")
set(BLAS "CBLAS" CACHE STRING "BLAS library to use")
set_property(CACHE BLAS PROPERTY STRINGS "CBLAS" "MKL")
set(PARPE_ENABLE_MPI TRUE CACHE BOOL "Use MPI?")
set(PARPE_ENABLE_FIDES FALSE CACHE BOOL "Enable fides optimizer?")
set(PARPE_ENABLE_IPOPT TRUE CACHE BOOL "Enable ipopt optimizer?")
set(PARPE_ENABLE_CERES TRUE CACHE BOOL "Enable ceres optimizer?")
set(PARPE_ENABLE_DLIB FALSE CACHE BOOL "Enable dlib optimizers?")
set(PARPE_ENABLE_TOMS611 FALSE CACHE BOOL "Enable toms611 optimizer?")
set(PARPE_ENABLE_FSQP FALSE CACHE BOOL "Enable FSQP optimizer?")
set(PARPE_EXPORT_PACKAGE FALSE CACHE BOOL
    "Export this build to CMake registry?")
set(PARPE_BUILD_OPTIMIZED TRUE CACHE BOOL "Optimize for current CPU?")

if(${PARPE_BUILD_OPTIMIZED})
    include(BuildOptimized)
endif(${PARPE_BUILD_OPTIMIZED})


# OpenMP?
find_package(OpenMP)

# Enable MPI?
if(${PARPE_ENABLE_MPI})
    find_package(MPI REQUIRED)
endif(${PARPE_ENABLE_MPI})
# </Build options>

if(${PARPE_ENABLE_FIDES})
    find_package(Fides REQUIRED)
endif(${PARPE_ENABLE_FIDES})

find_package(PkgConfig)
pkg_search_module(IPOPT REQUIRED IMPORTED_TARGET GLOBAL ipopt>=3.11.0 )

if(${PARPE_ENABLE_DLIB})
    set(DLIB_DIR "${CMAKE_CURRENT_LIST_DIR}/ThirdParty/dlib-19.7"
        CACHE PATH "DLIB base directory")
endif(${PARPE_ENABLE_DLIB})

# HDF5
#set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 COMPONENTS CXX C HL REQUIRED)
set(HDF5_LIBRARIES
    ${HDF5_HL_LIBRARIES} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})

# Guidelines Support Library (GSL-lite) directory
set(GSL_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/ThirdParty/gsl")
# Throw on contract violation instead of terminate to get more meaningful
#  message
add_definitions(-Dgsl_CONFIG_CONTRACT_VIOLATION_THROWS)

# Build tests?
if(${BUILD_TESTING})
    enable_testing()

    # --allow-run-as-root is required for running on shippable systems
    set(TESTS_MPIEXEC_COMMAND
        mpiexec --allow-run-as-root --oversubscribe -n 4 CACHE STRING "")

    # <Coverage>
    # Create test coverage reports?
    set(GCOV_REPORT FALSE CACHE BOOL "Create GCOV report?")
    set(GCOVR_REPORT FALSE CACHE BOOL "Create GCOVR report?")
    set(COVERAGE_DEPENDENCIES
        unittests_common unittests_amici unittests_loadbalancer
        unittests_optimization)

    if(${GCOV_REPORT})
        if(${CMAKE_COMPILER_IS_GNUCXX})
            # NOTE: LCOV uses absolute paths for exclusion
            #  arguments must be enclosed in single quotes
            set(COVERAGE_EXCLUDES
                '/usr/*' '${CMAKE_CURRENT_LIST_DIR}/ThirdParty/*'
                '${CMAKE_CURRENT_LIST_DIR}/deps/*'
                '${CMAKE_CURRENT_LIST_DIR}/examples/*')
            set(GCOV_LIBRARY gcov)
            include(CodeCoverage)
            APPEND_COVERAGE_COMPILER_FLAGS()
            setup_target_for_coverage(
                NAME ${PROJECT_NAME}_coverage
                EXECUTABLE ctest
                DEPENDENCIES ${COVERAGE_DEPENDENCIES})
        else()
            message(WARNING "Not creating coverage report, because compiler "
                "is not GCC but ${CMAKE_CXX_COMPILER_ID}.")
        endif()
    endif()

    if(${GCOVR_REPORT})
        if(${CMAKE_COMPILER_IS_GNUCXX})
            # NOTE: GCOVR uses regexs on relative paths from
            #  ${PROJECT_BINARY_DIR} for exclusion
            set(COVERAGE_EXCLUDES "^ThirdParty/" "^deps/" "^examples/")
            set(GCOV_LIBRARY gcov)
            include(CodeCoverage)
            APPEND_COVERAGE_COMPILER_FLAGS()
            setup_target_for_coverage_cobertura(
                NAME ${PROJECT_NAME}_coverage_cobertura
                EXECUTABLE ctest
                DEPENDENCIES ${COVERAGE_DEPENDENCIES})
        else()
            message(WARNING "Not creating coverage report, because compiler "
                "is not GCC but ${CMAKE_CXX_COMPILER_ID}.")
        endif()
    endif()
    # </Coverage>
endif(${BUILD_TESTING})


# <Documentation (doxygen)>
if(BUILD_DOCS)
    find_package(Doxygen)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(doxyfile_in ${CMAKE_CURRENT_LIST_DIR}/doc/Doxyfile.in)
    set(doxyfile ${CMAKE_CURRENT_LIST_DIR}/doc/Doxyfile)
    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
        COMMAND cd "${CMAKE_CURRENT_LIST_DIR}/doc/doxy/latex/"
            && pdflatex refman.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()
# </Documentation (doxygen)>

include(clang-tools)

# calc lines of code
add_custom_target(cloc
    COMMAND cloc misc src include examples --fullpath
        --not-match-d=amici/examples/steadystate/model/*
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    COMMENT "Counting lines of code"
    VERBATIM)

add_subdirectory(templates)

add_subdirectory(src)

if(${BUILD_TESTING})
    add_subdirectory(tests)
endif()

if(${BUILD_EXAMPLES})
    add_subdirectory(examples)
endif()

if(BUILD_PYTHON_MODULE AND CMAKE_BUILD_TYPE MATCHES Debug)
    # Only build swig/python module when we are in debug mode
    # This requires PIC, but we don't want PIC for release build at the moment
    add_subdirectory(swig)
endif()


# create single library
add_library(${PROJECT_NAME} /dev/null)
target_link_libraries(${PROJECT_NAME}
    PUBLIC parpecommon parpeloadbalancer parpeoptimization parpeamici)

# <Export cmake configuration>
install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
    )
export(EXPORT ParPETargets
    #DESTINATION "lib/cmake"
    NAMESPACE "Upstream::"
    FILE "ParPETargets.cmake"
    )
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ParPETargets.cmake"
    DESTINATION "lib/cmake")

set(LIB_INSTALL_DIR lib)
include(CMakePackageConfigHelpers)
configure_package_config_file(
    "cmake/ParPEConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/ParPEConfig.cmake"
    INSTALL_DESTINATION "${LIB_INSTALL_DIR}/cmake/")

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/ParPEConfigVersion.cmake"
    VERSION "${major}.${minor}.${patch}"
    COMPATIBILITY ExactVersion)

install(
    FILES
    "${CMAKE_CURRENT_BINARY_DIR}/ParPEConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/ParPEConfigVersion.cmake"
    DESTINATION "${LIB_INSTALL_DIR}/cmake/")


# Register package
if(${PARPE_EXPORT_PACKAGE})
    export(PACKAGE ParPE)
endif(${PARPE_EXPORT_PACKAGE})
# </Export cmake configuration>

# <Package configuration>
set (CPACK_PACKAGE_VERSION_MAJOR "${major}")
set (CPACK_PACKAGE_VERSION_MINOR "${minor}")
set (CPACK_PACKAGE_VERSION_PATCH "${patch}")
include (CPack)
# </Package configuration>
