set(AMICI_DIR "/home/dweindl/src/AMICI.github" CACHE PATH "AMICI root directory")

include_directories(
    ${AMICI_DIR}/include
    ../loadbalancer/src/
    ../common/src
    ../optimization/src
)

add_subdirectory(examples)
