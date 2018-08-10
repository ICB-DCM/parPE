find_package(Amici HINTS ${CMAKE_SOURCE_DIR}/deps/AMICI/build)

if(CMAKE_BUILD_TYPE MATCHES Release)
    set(Boost_USE_STATIC_LIBS TRUE)
else()
    # For python module we need -fPIC which is only the case with shared libs
    set(Boost_USE_STATIC_LIBS FALSE)
endif()
find_package(Boost COMPONENTS serialization REQUIRED)

add_subdirectory(src)

if(${BUILD_EXAMPLES})
    add_subdirectory(examples)
endif(${BUILD_EXAMPLES})

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif(${BUILD_TESTS})
