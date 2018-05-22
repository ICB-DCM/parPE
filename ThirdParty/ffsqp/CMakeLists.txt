project(fsqp)

add_library(${PROJECT_NAME} ffsqp.c ql0001.c)

target_link_libraries(${PROJECT_NAME}
    ${CMAKE_SOURCE_DIR}/ThirdParty/libf2c/libf2c.a
)

install(TARGETS ${PROJECT_NAME} EXPORT ParPETargets ARCHIVE DESTINATION lib)
