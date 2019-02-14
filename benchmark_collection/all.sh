#!/usr/bin/env bash
# Import, build, and run all benchmark models

BENCHMARK_DIR=${HOME}/src/Benchmark-Models-Caro/hackathon_contributions_new_data_format/

for MODEL in $(ls -1d ${BENCHMARK_DIR}/*/); do
    printf '=%.0s' {1..20}
    printf ${MODEL}
    printf '=%.0s' {1..20}
    echo
    ./import_and_run.sh ${MODEL}
    grep $(basename ${MODEL}) nllh.txt
    printf '=%.0s' {1..100}
    echo
done
