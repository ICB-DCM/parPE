set(AMICI_DIR "/home/dweindl/src/AMICI.github" CACHE PATH "AMICI root directory")

find_package(Boost COMPONENTS serialization)

if(APPLE)
    set(SHARED_OBJECT_EXTENSION .dylib)
else()
    set(SHARED_OBJECT_EXTENSION .so)
endif()

set(SUITESPARSE_DIR "${AMICI_DIR}/SuiteSparse/" CACHE PATH "SuiteSparse root directory")
set(SUITESPARSE_INCLUDE_DIRS "${AMICI_DIR}/SuiteSparse/include" CACHE PATH "SuiteSparse root directory")
set(SUITESPARSE_LIB_DIR "${AMICI_DIR}/SuiteSparse/lib" CACHE PATH "SuiteSparse library directory")
set(SUITESPARSE_LIBRARIES
    "${SUITESPARSE_LIB_DIR}/libcolamd${SHARED_OBJECT_EXTENSION}"
    "${SUITESPARSE_LIB_DIR}/libklu${SHARED_OBJECT_EXTENSION}"
    "${SUITESPARSE_LIB_DIR}/libbtf${SHARED_OBJECT_EXTENSION}"
    "${SUITESPARSE_LIB_DIR}/libamd${SHARED_OBJECT_EXTENSION}"
    "${SUITESPARSE_LIB_DIR}/libsuitesparseconfig${SHARED_OBJECT_EXTENSION}"
)

set(SUNDIALS_LIB_DIR "${AMICI_DIR}/sundials/build/lib" CACHE PATH "Sundials library directory")
set(SUNDIALS_INCLUDE_DIRS "${AMICI_DIR}/sundials/include" CACHE PATH "Sundials include directory")
set(SUNDIALS_LIBRARIES "${SUNDIALS_LIB_DIR}/libsundials_nvecserial.so"
                       "${SUNDIALS_LIB_DIR}/libsundials_cvodes.so"
)

include_directories(
    src
    ${AMICI_DIR}/include
    ${AMICI_DIR}
    ../loadbalancer/src/
    ../common/src
    ../optimization/src
    ${SUITESPARSE_INCLUDE_DIRS}
    ${SUNDIALS_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIR}
)

add_subdirectory(examples)
add_subdirectory(src)

if(${BUILD_TESTS})
    add_subdirectory(tests)
endif(${BUILD_TESTS})
