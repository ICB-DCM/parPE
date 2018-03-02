# Version number from git
set(PARPE_VERSION "unknown")
include("${CMAKE_CURRENT_LIST_DIR}/getVersion.cmake")
configure_file(${INFILE} ${OUTFILE} @ONLY)
message(STATUS "Buildung version ${PARPE_VERSION}")
