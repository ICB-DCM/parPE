find_package(Amici HINTS ${CMAKE_SOURCE_DIR}/deps/AMICI/build)

set(Boost_USE_STATIC_LIBS TRUE)
find_package(Boost COMPONENTS serialization REQUIRED)

add_subdirectory(src)

if(${BUILD_EXAMPLES})
    add_subdirectory(examples)
endif(${BUILD_EXAMPLES})

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif(${BUILD_TESTS})
