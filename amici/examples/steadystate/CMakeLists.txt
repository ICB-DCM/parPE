if(APPLE)
    set(SHARED_OBJECT_EXTENSION .dylib)
else()
    set(SHARED_OBJECT_EXTENSION .so)
endif()

set(MODEL_LIBRARIES "model_steadystate")
set(MODEL_INCLUDE_DIRS "model")

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
    ${MODEL_INCLUDE_DIRS}
    ${AMICI_DIR}
    ${SUITESPARSE_INCLUDE_DIRS}
    ${SUNDIALS_INCLUDE_DIRS}
)


# AMICI-generated model files
add_subdirectory(model)

# example_steadystate executable
project(example_steadystate)

set(SRC_LIST_EXE
    main.cpp
    steadystateProblem.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE})

target_link_libraries(${PROJECT_NAME}
    parpeoptimization
    parpeloadbalancer
    ${MODEL_LIBRARIES}
    ${IPOPT_LIBRARIES}
    ${SUNDIALS_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${CERES_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)
# /example_steadystate executable

# example_steadystate_parallel executable
project(example_steadystate_parallel)

set(SRC_LIST_EXE_PARALLEL
    main_parallel.cpp
    steadystateProblem.cpp
    steadystateProblemParallel.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE_PARALLEL})

target_link_libraries(${PROJECT_NAME}
    parpeoptimization
    parpeloadbalancer
    ${MODEL_LIBRARIES}
    ${IPOPT_LIBRARIES}
    ${SUNDIALS_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${CERES_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)
# /example_steadystate_parallel executable


# example_steadystate_multi executable

project(example_steadystate_multi)

set(SRC_LIST_EXE_PARALLEL
    main_multicondition.cpp
    SteadyStateMultiConditionProblem.cpp
    steadystateProblem.cpp
    steadystateProblemParallel.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST_EXE_PARALLEL})

target_link_libraries(${PROJECT_NAME}
    parpeoptimization
    parpeloadbalancer
    parpeamici
    ${MODEL_LIBRARIES}
    ${IPOPT_LIBRARIES}
    ${SUNDIALS_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_C_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${CERES_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)
# /example_steadystate_multi executable
