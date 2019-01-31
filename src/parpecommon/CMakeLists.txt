project(parpecommon)

# write config header
configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/parpeConfig.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/parpeConfig.h
    @ONLY)

set(SRC_LIST
    misc.cpp
    logging.cpp
    hdf5Misc.cpp
    parpeException.cpp
    model.cpp
    costFunction.cpp
    functions.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_include_directories(${PROJECT_NAME}
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../
    PUBLIC ${CMAKE_SOURCE_DIR}/include/
    PUBLIC ${HDF5_INCLUDE_DIRS}
    PUBLIC ${GSL_INCLUDE_DIRS}
    )

target_link_libraries(${PROJECT_NAME}
    PUBLIC ${HDF5_LIBRARIES}
)

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)

install(
    DIRECTORY
    ${CMAKE_SOURCE_DIR}/include/parpecommon
    DESTINATION
    include)
