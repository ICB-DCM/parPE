include_directories(
    src/
    ../loadbalancer/src
    tests/
    ${IPOPT_INCLUDE_DIRS}
    ${CERES_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

set(PARPE_CERES_MINIGLOG_REDIRECT FALSE CACHE BOOL "Custom treatment of Ceres output? Requires miniglog headers in ceres directory.")
if(${PARPE_CERES_MINIGLOG_REDIRECT})
    add_definitions(-DPARPE_CERES_MINIGLOG_REDIRECT)
endif(${PARPE_CERES_MINIGLOG_REDIRECT})

if(${PARPE_ENABLE_TOMS611})
    add_subdirectory(../ThirdParty/toms611/ toms611)
    include_directories(../ThirdParty/toms611/)
endif(${PARPE_ENABLE_TOMS611})

if(${PARPE_ENABLE_DLIB})
    find_package(dlib REQUIRED
        HINTS "${CMAKE_CURRENT_LIST_DIR}/../ThirdParty/dlib-19.7/build/install/lib/cmake/dlib/")
    include_directories(${dlib_INCLUDE_DIRS})
endif(${PARPE_ENABLE_DLIB})

if(${PARPE_ENABLE_FSQP})
    add_subdirectory(../ThirdParty/ffsqp ffsqp)
    include_directories(../ThirdParty/ffsqp/)
endif(${PARPE_ENABLE_FSQP})

add_subdirectory(src)

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif(${BUILD_TESTS})
