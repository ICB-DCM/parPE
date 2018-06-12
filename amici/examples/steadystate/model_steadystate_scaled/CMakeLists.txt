project(model_steadystate_scaled)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckCXXCompilerFlag)
set(MY_CXX_FLAGS -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable)
foreach(FLAG ${MY_CXX_FLAGS})
    unset(CUR_FLAG_SUPPORTED CACHE)
    CHECK_CXX_COMPILER_FLAG(${FLAG} CUR_FLAG_SUPPORTED)
    if(${CUR_FLAG_SUPPORTED})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
    endif()
endforeach(FLAG)

find_package(Amici HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build)

set(MODEL_DIR ${CMAKE_CURRENT_LIST_DIR})

set(SRC_LIST_LIB model_steadystate_scaled_J.cpp 
model_steadystate_scaled_JB.cpp 
model_steadystate_scaled_JDiag.cpp 
model_steadystate_scaled_JSparse.cpp 
model_steadystate_scaled_JSparseB.cpp 
model_steadystate_scaled_Jv.cpp 
model_steadystate_scaled_JvB.cpp 
model_steadystate_scaled_Jy.cpp 
model_steadystate_scaled_dJydsigma.cpp 
model_steadystate_scaled_dJydy.cpp 
model_steadystate_scaled_dwdp.cpp 
model_steadystate_scaled_dwdx.cpp 
model_steadystate_scaled_dxdotdp.cpp 
model_steadystate_scaled_dydx.cpp 
model_steadystate_scaled_dydp.cpp 
model_steadystate_scaled_dsigmaydp.cpp 
model_steadystate_scaled_qBdot.cpp 
model_steadystate_scaled_sigma_y.cpp 
model_steadystate_scaled_sxdot.cpp 
model_steadystate_scaled_w.cpp 
model_steadystate_scaled_x0.cpp 
model_steadystate_scaled_sx0.cpp 
model_steadystate_scaled_xBdot.cpp 
model_steadystate_scaled_xdot.cpp 
model_steadystate_scaled_y.cpp 
${MODEL_DIR}/wrapfunctions.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST_LIB})
add_library(model ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
    
target_link_libraries(${PROJECT_NAME}
    PUBLIC Upstream::amici
)

set(SRC_LIST_EXE main.cpp)

add_executable(simulate_${PROJECT_NAME} ${SRC_LIST_EXE})
    
target_link_libraries(simulate_${PROJECT_NAME} ${PROJECT_NAME})

if($ENV{ENABLE_GCOV_COVERAGE})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()

## SWIG
option(ENABLE_SWIG "Build swig/python library?" ON)
if(ENABLE_SWIG)
    add_subdirectory(swig)
endif()


# <Export cmake configuration>
include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Targets
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
export(EXPORT ${PROJECT_NAME}Targets FILE ${PROJECT_NAME}Config.cmake
    NAMESPACE Upstream::
    )
# </Export cmake configuration>

