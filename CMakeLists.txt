project(parPE)

cmake_minimum_required(VERSION 3.1.0)

if(POLICY CMP0056) # for CMake 3.2 and above
   cmake_policy (SET CMP0056 OLD)
endif()

set(CMAKE_DEBUG_POSTFIX "-dbg")

set(USE_MPI TRUE CACHE BOOL "Use MPI?")
set(BUILD_TESTS TRUE CACHE BOOL "Build tests?")

set(BLAS "CBLAS" CACHE STRING "BLAS library to use")
set_property(CACHE BLAS PROPERTY STRINGS "CBLAS" "MKL")

if(${BLAS} STREQUAL "MKL")
    add_definitions(-DAMICI_BLAS_MKL)
    set(BLAS_LIBRARIES -lmkl CACHE STRING "")
    set(BLAS_INCLUDE_DIRS "" CACHE STRING "")
else()
    add_definitions(-DAMICI_BLAS_CBLAS)
    set(BLAS_INCLUDE_DIRS "" CACHE STRING "")
    set(BLAS_LIBRARIES -lcblas CACHE STRING "")
endif()

if(${USE_MPI})
    find_package(MPI REQUIRED)
    add_definitions(-DUSE_MPI)

    #FIXME
    #execute_process(COMMAND "${MPI_C_COMPILER} --showme:incdirs" OUTPUT_VARIABLE MPI_INCLUDE_DIR)
    set(MPI_INCLUDE_DIR "/usr/lib/openmpi/include/" CACHE PATH "MPI include directory")
endif(${USE_MPI})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wno-unused-function -fopenmp -D_GNU_SOURCE") # -D_GNU_SOURCE for pthread recursive mutex issues
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-unused-function -fopenmp -D_GNU_SOURCE") # -D_GNU_SOURCE for pthread recursive mutex issues
# -std=c++11 or 14 break compilation due to issues with hdf5/mpi/stl::map

execute_process(COMMAND sh -c "cd ${CMAKE_CURRENT_SOURCE_DIR} && git describe --abbrev=4 --dirty=-dirty --always --tags  | tr -d '\n'" OUTPUT_VARIABLE GIT_VERSION)
message(STATUS "Building version ${GIT_VERSION}")
add_definitions(-DGIT_VERSION="${GIT_VERSION}")

add_definitions(-DAMICI_WITHOUT_MATLAB)

set(IPOPT_DIR "/home/dweindl/src/_libs/Ipopt-3.12.6/" CACHE PATH "IpOpt root directory")
set(IPOPT_INCLUDE_DIRS "${IPOPT_DIR}/include/coin/"  CACHE PATH "IpOpt include directory")
set(IPOPT_LIBRARIES "${IPOPT_DIR}/lib/libipopt.so"  CACHE FILEPATH "IpOpt library")

set(CERES_INCLUDE_DIRS "/home/dweindl/src/_libs/ceres-solver/build/install/include/" "/usr/include/eigen3" CACHE PATH "CERES include directories")
set(CERES_LIBRARIES "/home/dweindl/src/_libs/ceres-solver/build/install/lib/libceres.so" CACHE FILEPATH "CERES libraries")

set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 COMPONENTS C HL REQUIRED)

find_package(Threads REQUIRED)

include_directories(
    common/src/
    ${IPOPT_INCLUDE_DIRS}
    ${CERES_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${MPI_INCLUDE_DIR}
    ${BLAS_INCLUDE_DIRS}
)

if(${BUILD_TESTS})
    enable_testing ()

    set(CPPUTEST_DIR "/home/dweindl/src/cpputest/" CACHE PATH "")

    set(CPPUTEST_INCLUDE_DIR "${CPPUTEST_DIR}/include")
    set(CPPUTEST_LIBRARY "${CPPUTEST_DIR}/cpputest_build/lib/libCppUTest.a" "${CPPUTEST_DIR}/cpputest_build/lib/libCppUTestExt.a")
    include_directories(${CPPUTEST_INCLUDE_DIR})



endif(${BUILD_TESTS})

add_subdirectory(common)
add_subdirectory(optimization)
add_subdirectory(loadbalancer)
