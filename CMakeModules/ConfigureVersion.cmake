# Version number from git
set(PARPE_VERSION "unknown")
find_package(Git)
if(Git_FOUND)
    execute_process(COMMAND sh -c "${GIT_EXECUTABLE} describe --abbrev=4 --dirty=-dirty --always --tags  | tr -d '\n'"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE PARPE_VERSION
        )
endif()
configure_file(${INFILE} ${OUTFILE} @ONLY)
message(STATUS "Buildung version ${PARPE_VERSION}")
