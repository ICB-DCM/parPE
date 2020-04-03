set (CTEST_PROJECT_NAME "parPE")
set(MemoryCheckCommandOptions
    "${MemoryCheckCommandOptions}
        -v
        --leak-check=full
        --show-leak-kinds=all
        --error-exitcode=1
        --gen-suppressions=all
        --suppressions=${CMAKE_CURRENT_LIST_DIR}/tests/parpe.supp")
