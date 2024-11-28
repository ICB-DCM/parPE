# Version number from git
set(PARPE_VERSION "unknown")
include("${CMAKE_CURRENT_LIST_DIR}/getVersion.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/split_version.cmake")
split_version(${PARPE_VERSION} libname major minor patch)
# message(STATUS "${libname} ${major} ${minor} ${patch}  ")
configure_file(${INFILE} ${OUTFILE} @ONLY)
message(STATUS "Building version ${PARPE_VERSION}")
