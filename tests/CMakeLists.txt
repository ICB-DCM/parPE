add_subdirectory(parpecommon)

if(${PARPE_ENABLE_MPI})
    add_subdirectory(parpeloadbalancer)
endif()

add_subdirectory(parpeoptimization)
add_subdirectory(parpeamici)
