project(parpecommon)

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

target_link_libraries(${PROJECT_NAME}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

get_target_property(PARPECOMMON_INCLUDE_DIRS ${PROJECT_NAME} INCLUDE_DIRECTORIES)
set_property(GLOBAL PROPERTY PARPECOMMON_INCLUDE_DIRS "${PARPECOMMON_INCLUDE_DIRS}")

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)
