cmake_minimum_required(VERSION 3.6)
cmake_policy(VERSION 3.6)

project(parpe)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/CMakeModules/)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(BuildType) # Ensure CMAKE_BUILD_TYPE is always set
include(BuildOptimized)

set(CMAKE_DEBUG_POSTFIX "-dbg")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wno-unused-function -fopenmp -D_GNU_SOURCE") # -D_GNU_SOURCE for pthread recursive mutex issues
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unused-function -fopenmp -D_GNU_SOURCE") # -D_GNU_SOURCE for pthread recursive mutex issues

if(BUILD_PYTHON_MODULE)
    # Build PIC code to be used for swig/python module
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()


# <versioning>
# Version number from git
add_custom_target(get_version ALL
    COMMAND ${CMAKE_COMMAND} -DINFILE=${PROJECT_SOURCE_DIR}/src/parpecommon/parpeVersion.h.in
    -DOUTFILE=${CMAKE_BINARY_DIR}/src/parpecommon/parpeVersion.h
    -P ${PROJECT_SOURCE_DIR}/CMakeModules/ConfigureVersion.cmake)

include("getVersion")
include("split_version")
split_version(${PARPE_VERSION} libname major minor patch)
# </versioning>


# <Build options>
set(BUILD_TESTS FALSE CACHE BOOL "Build tests? (requires googletest)")
set(BUILD_EXAMPLES TRUE CACHE BOOL "Build examples?")
set(BUILD_DOCS FALSE CACHE BOOL "Generate documentation?")
set(BUILD_PYTHON_MODULE FALSE CACHE BOOL "Build Python module?")
set(BLAS "CBLAS" CACHE STRING "BLAS library to use")
set_property(CACHE BLAS PROPERTY STRINGS "CBLAS" "MKL")
set(PARPE_ENABLE_MPI TRUE CACHE BOOL "Use MPI?")
set(IPOPT_DIR "" CACHE PATH "IpOpt root directory")
set(IPOPT_INCLUDE_DIRS "${IPOPT_DIR}/include/coin/"  CACHE PATH "IpOpt include directory")
set(IPOPT_LIBRARIES ${IPOPT_DIR}/lib/libipopt.a ${IPOPT_DIR}/lib/libcoinhsl.a gfortran CACHE STRINGS "IpOpt library")
set(CERES_INCLUDE_DIRS "" "/usr/include/eigen3" CACHE PATH "CERES include directories")
set(CERES_LIBRARIES "" CACHE STRINGS "CERES libraries")
set(PARPE_ENABLE_IPOPT TRUE CACHE BOOL "Enable ipopt optimizer?")
set(PARPE_ENABLE_CERES TRUE CACHE BOOL "Enable ceres optimizer?")
set(PARPE_ENABLE_DLIB FALSE CACHE BOOL "Enable dlib optimizers?")
set(PARPE_ENABLE_TOMS611 FALSE CACHE BOOL "Enable toms611 optimizer?")
set(PARPE_ENABLE_FSQP FALSE CACHE BOOL "Enable FSQP optimizer?")
set(PARPE_EXPORT_PACKAGE FALSE CACHE BOOL "Export this build to CMake registry?")

# Enable MPI?
if(${PARPE_ENABLE_MPI})
    find_package(MPI REQUIRED)
    if(IS_DIRECTORY ${MPI_INCLUDE_DIRS})
        # may be undefined. will be provided by MPI compiler wrapped.
        # this is only added for convenient IDE integration
        include_directories(${MPI_INCLUDE_DIRS})
    endif()
endif(${PARPE_ENABLE_MPI})
# </Build options>

if(${ENABLE_DLIB})
    set(DLIB_DIR "${CMAKE_CURRENT_LIST_DIR}/ThirdParty/dlib-19.7" CACHE FILEPATH "DLIB base directory")
endif(${ENABLE_DLIB})

# PThreads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)


# HDF5
#set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 COMPONENTS CXX C HL REQUIRED)
set(HDF5_LIBRARIES ${HDF5_HL_LIBRARIES} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})

# Guidelines Support Library (GSL-lite) directory
set(GSL_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/ThirdParty/gsl)

# Build tests?
if(${BUILD_TESTS})
    enable_testing()
    include (CTest)

    # --allow-run-as-root is required for running on shippable systems
    set(TESTS_MPIEXEC_COMMAND mpiexec --allow-run-as-root -n 4 CACHE STRING "")

    # googletest
    add_definitions(-DGTEST_LANGUAGE_CXX11)
    find_package(Threads REQUIRED)
    if ($ENV{GOOGLETEST_DIR})
        set(GOOGLETEST_DIR $ENV{GOOGLETEST_DIR})
    else()
        set(GOOGLETEST_DIR ${CMAKE_CURRENT_LIST_DIR}/ThirdParty/googletest)
    endif()
    if (EXISTS ${GOOGLETEST_DIR})
        set(GTestSrc ${GOOGLETEST_DIR}/googletest)
        set(GMockSrc ${GOOGLETEST_DIR}/googlemock)
    else ()
        message(FATAL_ERROR "No googletest src dir found - set GOOGLETEST_DIR "
            "or run ../ThirdParty/installGoogleTest.sh to enable!"
            " (Or disable building tests with BUILD_TESTS=OFF)")
    endif ()
    include_directories(${GTestSrc} ${GTestSrc}/include ${GMockSrc} ${GMockSrc}/include)

    # <Coverage>
    # Create test coverage reports?
    set(GCOV_REPORT FALSE CACHE BOOL "Create GCOV report?")
    set(GCOVR_REPORT FALSE CACHE BOOL "Create GCOVR report?")
    set(COVERAGE_DEPENDENCIES unittests_common unittests_amici unittests_loadbalancer unittests_optimization)

    if(${GCOV_REPORT})
        if(${CMAKE_COMPILER_IS_GNUCXX})
            # NOTE: LCOV uses absolute paths for exclusion; arguments must be enclosed in single quotes
            set(COVERAGE_EXCLUDES '/usr/*' '${CMAKE_CURRENT_LIST_DIR}/ThirdParty/*' '${CMAKE_CURRENT_LIST_DIR}/deps/*' '${CMAKE_CURRENT_LIST_DIR}/examples/*')
            set(GCOV_LIBRARY gcov)
            include(CodeCoverage)
            APPEND_COVERAGE_COMPILER_FLAGS()
            setup_target_for_coverage(
                NAME ${PROJECT_NAME}_coverage
                EXECUTABLE ctest -R unittests
                DEPENDENCIES ${COVERAGE_DEPENDENCIES})
        else()
            message("Not creating coverage report, because compiler is not GCC but ${CMAKE_CXX_COMPILER_ID}.")
        endif()
    endif()

    if(${GCOVR_REPORT})
        if(${CMAKE_COMPILER_IS_GNUCXX})

            # NOTE: GCOVR uses regexs on relative paths from ${PROJECT_BINARY_DIR} for exclusion
            set(COVERAGE_EXCLUDES "^ThirdParty/" "^deps/" "^examples/")
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
    # </Coverage>
endif(${BUILD_TESTS})


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
        COMMAND cd ${CMAKE_CURRENT_LIST_DIR}/doc/doxy/latex/ && pdflatex refman.tex
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()
# </Documentation (doxygen)>

include(clang-tools)

# calc lines of code
add_custom_target(cloc
    COMMAND cloc misc src include examples --fullpath --not-match-d=amici/examples/steadystate/model/*
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    COMMENT "Counting lines of code"
    VERBATIM)

add_subdirectory(templates)

add_subdirectory(src)

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif()

if(${BUILD_EXAMPLES})
    add_subdirectory(examples)
endif()

if(BUILD_PYTHON_MODULE AND CMAKE_BUILD_TYPE MATCHES Debug)
    # Only build swig/python module when we are in debug mode
    # This requires PIC, but we don't want PIC for release build at the moment
    add_subdirectory(python/swig)
endif()


# create single library
add_library(${PROJECT_NAME} /dev/null)
target_link_libraries(${PROJECT_NAME} PUBLIC parpecommon parpeloadbalancer parpeoptimization parpeamici)


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
    CMakeModules/ParPEConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/ParPEConfig.cmake
    INSTALL_DESTINATION ${LIB_INSTALL_DIR}/cmake/)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/ParPEConfigVersion.cmake
    VERSION "${major}.${minor}.${patch}"
    COMPATIBILITY ExactVersion )

install(
    FILES
    ${CMAKE_CURRENT_BINARY_DIR}/ParPEConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/ParPEConfigVersion.cmake
    DESTINATION ${LIB_INSTALL_DIR}/cmake/)

configure_file(CMakeModules/ParPEConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/ParPEConfig.cmake"
    @ONLY
    )

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
