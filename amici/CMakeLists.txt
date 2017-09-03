find_package(Amici)

find_package(Boost COMPONENTS serialization)

get_property(AmiciConfigIncludes TARGET Upstream::amici PROPERTY INCLUDE_DIRECTORIES)

include_directories(
    src
    ../loadbalancer/src/
    ../common/src
    ../optimization/src
    ${Boost_INCLUDE_DIR}
    ${AmiciConfigIncludes}
)

add_subdirectory(src)

if(${BUILD_EXAMPLES})
    add_subdirectory(examples)
endif(${BUILD_EXAMPLES})

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif(${BUILD_TESTS})
