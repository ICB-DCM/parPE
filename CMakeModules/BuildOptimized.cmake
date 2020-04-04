# Enable all native optimizations in Release builds
if(CMAKE_BUILD_TYPE MATCHES Release OR CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG(-xHOST COMPILER_SUPPORTS_XHOST)
    CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCHNATIVE)

    if(COMPILER_SUPPORTS_XHOST)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST")
    endif()
    if(COMPILER_SUPPORTS_MARCHNATIVE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif()
endif()
