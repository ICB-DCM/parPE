set(CTEST_PROJECT_NAME "parPE")
set(MEMORYCHECK_COMMAND_OPTIONS
    "${MEMORYCHECK_COMMAND_OPTIONS} \
        -v \
        --leak-check=full \
        --show-leak-kinds=all \
        --error-exitcode=1 \
        --gen-suppressions=all \
        --suppressions=${CMAKE_CURRENT_LIST_DIR}/tests/parpe.supp")
