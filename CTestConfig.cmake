set(MEMORYCHECK_COMMAND_OPTIONS "${MEMORYCHECK_COMMAND_OPTIONS} --leak-check=full --error-exitcode=1 --gen-suppressions=all -v")
set(MEMORYCHECK_SUPPRESSIONS_FILE ${CMAKE_CURRENT_LIST_DIR}/unittests_common.supp)
