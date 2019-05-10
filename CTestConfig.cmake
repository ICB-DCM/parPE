set (CTEST_PROJECT_NAME "parPE")
set(MEMORYCHECK_COMMAND_OPTIONS "${MEMORYCHECK_COMMAND_OPTIONS} --leak-check=full --error-exitcode=1 --gen-suppressions=all -v --suppressions=${CMAKE_CURRENT_LIST_DIR}/tests/parpe.supp")
