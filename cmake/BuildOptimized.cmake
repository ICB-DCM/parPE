# Enable all native optimizations in Release builds
if(CMAKE_BUILD_TYPE MATCHES Release OR CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-xHOST COMPILER_SUPPORTS_XHOST)
  check_cxx_compiler_flag(-march=native COMPILER_SUPPORTS_MARCHNATIVE)

  if(COMPILER_SUPPORTS_XHOST)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHOST")
  elseif(COMPILER_SUPPORTS_MARCHNATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif()
endif()
