project(toms611)

add_library(${PROJECT_NAME} toms611.cpp)

target_link_libraries(${PROJECT_NAME}
    ${CMAKE_CURRENT_LIST_DIR}/../libf2c/libf2c.a
)

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets ARCHIVE DESTINATION lib)
