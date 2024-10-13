project(parpeloadbalancer)

set(SRC_LIST loadBalancerWorker.cpp loadBalancerMaster.cpp)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(${PROJECT_NAME} PUBLIC parpecommon)

if(${PARPE_ENABLE_MPI})
  target_compile_options(MPI::MPI_C INTERFACE -Wno-unused-parameter -Wno-error)

  target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC MPI::MPI_CXX
    PUBLIC MPI::MPI_C)
endif(${PARPE_ENABLE_MPI})

install(
  TARGETS ${PROJECT_NAME}
  EXPORT ParPETargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/parpeloadbalancer
        DESTINATION include)
