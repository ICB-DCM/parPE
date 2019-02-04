project(unittests_common)

add_test (NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} -c)

set(SRC_LIST
    main.cpp
    commonTests.cpp
    hdf5MiscTests.cpp
    ${CMAKE_SOURCE_DIR}/tests/parpecommon/testingMisc.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST})

#target_compile_options(${PROJECT_NAME} PUBLIC
#    -include ${CMAKE_CURRENT_LIST_DIR}/../src/STLCompatibleMemoryLeakDetectorMacros.h)

target_link_libraries(${PROJECT_NAME}
    ${CppUTest_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    parpecommon
    ${GCOV_LIBRARY}
)
