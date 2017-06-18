project(parpecommon)

set(SRC_LIST
    misc.c
    logging.c
    hdf5Misc.c
)

add_library(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
)
