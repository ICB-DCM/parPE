project(parpe)

# Stick to version 2.8.12 for now to work on ICB servers
cmake_minimum_required(VERSION 2.8.12)

if(POLICY CMP0056) # for CMake 3.2 and above
    cmake_policy (SET CMP0056 OLD)
endif()
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/CMakeModules/)
set(CMAKE_DEBUG_POSTFIX "-dbg")
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wno-unused-function -fopenmp -D_GNU_SOURCE") # -D_GNU_SOURCE for pthread recursive mutex issues
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unused-function -fopenmp -D_GNU_SOURCE") # -D_GNU_SOURCE for pthread recursive mutex issues
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

include("getVersion")
include("split_version")
split_version(${PARPE_VERSION} libname major minor patch)
message(STATUS "${libname} ${major} ${minor} ${patch}  ")

# Build options
set(BUILD_TESTS TRUE CACHE BOOL "Build tests? (requires cppUTest)")
set(BUILD_EXAMPLES TRUE CACHE BOOL "Build examples?")
set(BUILD_DOCS FALSE CACHE BOOL "Generate documentation?")
set(BLAS "CBLAS" CACHE STRING "BLAS library to use")
set_property(CACHE BLAS PROPERTY STRINGS "CBLAS" "MKL")
#set(USE_MPI TRUE CACHE BOOL "Use MPI?")
set(IPOPT_DIR "" CACHE PATH "IpOpt root directory")
set(IPOPT_INCLUDE_DIRS "${IPOPT_DIR}/include/coin/"  CACHE PATH "IpOpt include directory")
set(IPOPT_LIBRARIES ${IPOPT_DIR}/lib/libipopt.a ${IPOPT_DIR}/lib/libcoinhsl.a gfortran CACHE STRINGS "IpOpt library")
set(CERES_INCLUDE_DIRS "" "/usr/include/eigen3" CACHE PATH "CERES include directories")
set(CERES_LIBRARIES "" CACHE STRINGS "CERES libraries")
set(ENABLE_DLIB FALSE CACHE BOOL "Enable dlib optimizers?")
set(ENABLE_TOMS611 FALSE CACHE BOOL "Enable toms611 optimizers?")
set(ENABLE_FSQP FALSE CACHE BOOL "Enable FSQP optimizers?")

# Enable MPI? Currently mandatory.
set(USE_MPI ON)
if(${USE_MPI})
    find_package(MPI REQUIRED)
    add_definitions(-DUSE_MPI)
    set(MPI_INCLUDE_DIR "/usr/lib/openmpi/include/" CACHE PATH "MPI include directory")
endif(${USE_MPI})

if(${ENABLE_DLIB})
    set(DLIB_DIR "${CMAKE_CURRENT_LIST_DIR}/ThirdParty/dlib-19.7" CACHE FILEPATH "DLIB base directory")
    add_definitions(-DPARPE_DLIB_ENABLED)
endif(${ENABLE_DLIB})

if(${ENABLE_TOMS611})
    add_definitions(-DPARPE_TOMS611_ENABLED)
endif(${ENABLE_TOMS611})

if(${ENABLE_FSQP})
    add_definitions(-DPARPE_FSQP_ENABLED)
endif(${ENABLE_FSQP})

# PThreads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)


# Version number from git
add_custom_target(get_version ALL
    COMMAND ${CMAKE_COMMAND} -DINFILE=${PROJECT_SOURCE_DIR}/common/src/parpeVersion.h.in
    -DOUTFILE=${PROJECT_SOURCE_DIR}/common/src/parpeVersion.h
    -P ${PROJECT_SOURCE_DIR}/CMakeModules/ConfigureVersion.cmake)

# HDF5
#set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 COMPONENTS CXX C HL REQUIRED)


include_directories(
    common/src/
    ${HDF5_INCLUDE_DIRS}
    ${MPI_INCLUDE_DIR}
    )


# Build tests?
if(${BUILD_TESTS})
    include (CTest)

    # --allow-run-as-root is required for running on shippable systems
    set(TESTS_MPIEXEC_COMMAND mpiexec --allow-run-as-root -n 4 CACHE STRING "")

    # CppUTest
    add_definitions(-DD_MemoryLeakWarningPlugin_h) # prevent undefined reference to operator new; new is redefined even when memoryleakdetector is disabled
    set(CPPUTEST_DIR "${CMAKE_CURRENT_LIST_DIR}/deps/AMICI/ThirdParty/cpputest-master/build-noleakcheck/install/" CACHE PATH "cppUTest base directory")
    set(CPPUTEST_INCLUDE_DIR "${CPPUTEST_DIR}/include")
    set(CPPUTEST_LIBRARY "${CPPUTEST_DIR}/lib/libCppUTest.a" "${CPPUTEST_DIR}/lib/libCppUTestExt.a")
    include_directories(${CPPUTEST_INCLUDE_DIR})

    # Create test coverage reports?
    set(GCOV_REPORT FALSE CACHE BOOL "Create GCOV report?")
    set(GCOVR_REPORT FALSE CACHE BOOL "Create GCOVR report?")
    set(COVERAGE_DEPENDENCIES unittests_common unittests_amici unittests_loadbalancer unittests_optimization)

    if(${GCOV_REPORT})
        if(${CMAKE_COMPILER_IS_GNUCXX})
            # NOTE: LCOV uses absolute paths for exclusion; arguments must be enclosed in single quotes
            set(COVERAGE_EXCLUDES '/usr/*' '${CMAKE_CURRENT_LIST_DIR}/ThirdParty/*' '${CMAKE_CURRENT_LIST_DIR}/deps/*' '${CMAKE_CURRENT_LIST_DIR}/amici/examples/steadystate/model/*')
            set(GCOV_LIBRARY gcov)
            include(CodeCoverage)
            APPEND_COVERAGE_COMPILER_FLAGS()
            setup_target_for_coverage(
                NAME ${PROJECT_NAME}_coverage
                EXECUTABLE ctest
                DEPENDENCIES ${COVERAGE_DEPENDENCIES})
        else()
            message("Not creating coverage report, because coamici/examples/steadystate/modelamici/examples/steadystate/modelmpiler is not GCC but ${CMAKE_CXX_COMPILER_ID}.")
        endif()
    endif()

    if(${GCOVR_REPORT})
        if(${CMAKE_COMPILER_IS_GNUCXX})

            # NOTE: GCOVR uses regexs on relative paths from ${PROJECT_BINARY_DIR} for exclusion
            set(COVERAGE_EXCLUDES "^ThirdParty/" "^deps/" "^amici/examples/steadystate/model/")
            set(GCOV_LIBRARY gcov)
            include(CodeCoverage)
            APPEND_COVERAGE_COMPILER_FLAGS()
            setup_target_for_coverage_cobertura(
                NAME ${PROJECT_NAME}_coverage_cobertura
                EXECUTABLE ctest
                DEPENDENCIES ${COVERAGE_DEPENDENCIES})
        else()
            message("Not creating coverage report, because compiler is not GCC but ${CMAKE_CXX_COMPILER_ID}.")
        endif()
    endif()

endif(${BUILD_TESTS})


# Documentation (doxygen)
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
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()

include(CMakeModules/clang-tools.cmake)

add_custom_target(cloc
    COMMAND cloc misc amici common optimization loadbalancer --fullpath --not-match-d=amici/examples/steadystate/model/*
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    COMMENT "Counting lines of code"
    VERBATIM)


add_subdirectory(common)
add_subdirectory(optimization)
add_subdirectory(loadbalancer)
add_subdirectory(amici)

# create single library
add_library(${PROJECT_NAME} "")
target_link_libraries(${PROJECT_NAME} parpecommon parpeloadbalancer parpeoptimization parpeamici)

# <Export cmake configuration>
include(GenerateExportHeader)
generate_export_header(parpe)
install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
    )
export(EXPORT ParPETargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/ParPETargets.cmake"
    NAMESPACE Upstream::
    )
get_property(ParPEConfigIncludes DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
get_property(ParPEConfigIncludesCommon GLOBAL PROPERTY PARPECOMMON_INCLUDE_DIRS)
get_property(ParPEConfigIncludesLoadbalancer GLOBAL PROPERTY PARPELOADBALANCER_INCLUDE_DIRS)
get_property(ParPEConfigIncludesOptimization GLOBAL PROPERTY PARPEOPTIMIZATION_INCLUDE_DIRS)
get_property(ParPEConfigIncludesAmici GLOBAL PROPERTY PARPEAMICI_INCLUDE_DIRS)

configure_file(CMakeModules/ParPEConfig.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/ParPEConfig.cmake"
    @ONLY
    )
# Register package
export(PACKAGE ParPE)
# </Export cmake configuration>
