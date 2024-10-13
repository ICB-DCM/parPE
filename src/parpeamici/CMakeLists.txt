# dependency: AMICI
find_package(Amici PATHS ${CMAKE_SOURCE_DIR}/deps/AMICI/build)
if(Amici_FOUND)
  message(STATUS "Found AMICI version ${Amici_VERSION} at ${Amici_DIR}")
endif()

# dependency: BOOST
if(CMAKE_BUILD_TYPE MATCHES Release)
  set(Boost_USE_STATIC_LIBS TRUE)
else()
  # For python module we need -fPIC which is only the case with shared libs
  set(Boost_USE_STATIC_LIBS FALSE)
endif()
find_package(
  Boost
  COMPONENTS serialization thread
  REQUIRED)

project(parpeamici)

set(SRC_LIST
    multiConditionDataProvider.cpp
    multiConditionProblem.cpp
    optimizationApplication.cpp
    amiciSimulationRunner.cpp
    simulationResultWriter.cpp
    standaloneSimulator.cpp
    amiciMisc.cpp
    hierarchicalOptimization.cpp
    hierarchicalOptimizationAnalyticalParameterProvider.cpp)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_include_directories(
  ${PROJECT_NAME}
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC parpeloadbalancer
  PUBLIC parpecommon
  PUBLIC Upstream::amici
  PUBLIC Boost::serialization
  # For not fully understood reasons, `parpeoptimization` needs to come last.
  # There was an issue with the Ubuntu 24.04 ipopt package, which is pulling in
  # libdmumps_seq-5.6.so, which exports dummy MPI_* symbols that cause problems
  # when building parpe with MPI. Run-time resolution would pick up the
  # incorrect MPI symbols. We need to make sure that MPI::MPI_C linker flags
  # come first, but CMake first lists all direct dependencies, and only then the
  # transitive dependencies (in this case, in an unfortunate order).
  # INTERFACE_LINK_LIBRARIES_DIRECT did not help. Note that in this case, parpe
  # should *not* be compiled with the MPI compiler wrappers (mpicc, mpicxx),
  # because this will prevent CMake from adding the MPI linker flags.
  PUBLIC parpeoptimization)

if(${OpenMP_FOUND})
  target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_CXX)
endif()

install(
  TARGETS ${PROJECT_NAME}
  EXPORT ParPETargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/parpeamici DESTINATION include)
