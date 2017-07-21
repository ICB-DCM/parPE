project(parPE)

# Stick to version 2.8.12 for now to work on ICB servers
cmake_minimum_required(VERSION 2.8.12)

if(POLICY CMP0056) # for CMake 3.2 and above
   cmake_policy (SET CMP0056 OLD)
endif()

set(CMAKE_DEBUG_POSTFIX "-dbg")
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wno-unused-function -fopenmp -D_GNU_SOURCE") # -D_GNU_SOURCE for pthread recursive mutex issues
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unused-function -fopenmp -D_GNU_SOURCE") # -D_GNU_SOURCE for pthread recursive mutex issues
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)


# Build options
set(BUILD_TESTS TRUE CACHE BOOL "Build tests? (requires cppUTest)")
set(BUILD_EXAMPLES TRUE CACHE BOOL "Build examples?")
set(BUILD_DOCS FALSE CACHE BOOL "Generate documentation?")
set(BLAS "CBLAS" CACHE STRING "BLAS library to use")
set_property(CACHE BLAS PROPERTY STRINGS "CBLAS" "MKL")
set(USE_MPI TRUE CACHE BOOL "Use MPI?")
set(IPOPT_DIR "" CACHE PATH "IpOpt root directory")
set(IPOPT_INCLUDE_DIRS "${IPOPT_DIR}/include/coin/"  CACHE PATH "IpOpt include directory")
set(IPOPT_LIBRARIES "${IPOPT_DIR}/lib/libipopt.so"  CACHE FILEPATH "IpOpt library")
set(CERES_INCLUDE_DIRS "" "/usr/include/eigen3" CACHE PATH "CERES include directories")
set(CERES_LIBRARIES "" CACHE FILEPATH "CERES libraries")

# Which BLAS to use?
if(${BLAS} STREQUAL "MKL")
    add_definitions(-DAMICI_BLAS_MKL)
    set(BLAS_LIBRARIES -lmkl CACHE STRING "")
    set(BLAS_INCLUDE_DIRS "" CACHE STRING "")
else()
    add_definitions(-DAMICI_BLAS_CBLAS)
    set(BLAS_INCLUDE_DIRS "" CACHE STRING "")
    set(BLAS_LIBRARIES -lcblas CACHE STRING "")
endif()


# Enable MPI? Currently mandatory.
if(${USE_MPI})
    find_package(MPI REQUIRED)
    add_definitions(-DUSE_MPI)
    set(MPI_INCLUDE_DIR "/usr/lib/openmpi/include/" CACHE PATH "MPI include directory")
endif(${USE_MPI})

# PThreads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)


# Version number from git
execute_process(COMMAND sh -c "cd ${CMAKE_CURRENT_SOURCE_DIR} && git describe --abbrev=4 --dirty=-dirty --always --tags  | tr -d '\n'" OUTPUT_VARIABLE GIT_VERSION)
message(STATUS "Building version ${GIT_VERSION}")
add_definitions(-DGIT_VERSION="${GIT_VERSION}")


# HDF5
set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 COMPONENTS C HL REQUIRED)


include_directories(
    common/src/
    ${IPOPT_INCLUDE_DIRS}
    ${CERES_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${MPI_INCLUDE_DIR}
    ${BLAS_INCLUDE_DIRS}
)


# Build tests?
if(${BUILD_TESTS})
    enable_testing ()

    set(CPPUTEST_DIR "" CACHE PATH "cppUTest base directory")
    set(CPPUTEST_INCLUDE_DIR "${CPPUTEST_DIR}/include")
    set(CPPUTEST_LIBRARY "${CPPUTEST_DIR}/cpputest_build/lib/libCppUTest.a" "${CPPUTEST_DIR}/cpputest_build/lib/libCppUTestExt.a")
    include_directories(${CPPUTEST_INCLUDE_DIR})

    # Create test coverage reports?
    set(GCOV_REPORT FALSE CACHE BOOL "Create GCOV report?")
    set(GCOVR_REPORT FALSE CACHE BOOL "Create GCOVR report?")

    if(${CMAKE_COMPILER_IS_GNUCXX})
        set(COVERAGE_EXCLUDES 'boost/*' '/usr/include/*' 'cxx/*', 'ThirdParty/*')
        if(${GCOV_REPORT})
            set(GCOV_LIBRARY gcov)
            include(CodeCoverage)
            APPEND_COVERAGE_COMPILER_FLAGS()
            setup_target_for_coverage(
                NAME ${PROJECT_NAME}_coverage
                EXECUTABLE ctest
                DEPENDENCIES parpeoptimization)
        endif()
        if(${GCOVR_REPORT})
            set(GCOV_LIBRARY gcov)
            include(CodeCoverage)
            APPEND_COVERAGE_COMPILER_FLAGS()
            setup_target_for_coverage_cobertura(
                NAME ${PROJECT_NAME}_coverage_cobertura
                EXECUTABLE ctest
                DEPENDENCIES parpeoptimization)
        endif()
    else()
        message("Not creating coverage report, because compiler is not GCC but ${CMAKE_CXX_COMPILER_ID}.")
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


add_subdirectory(common)
add_subdirectory(optimization)
add_subdirectory(loadbalancer)
add_subdirectory(amici)
