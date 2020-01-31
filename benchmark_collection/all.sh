#!/usr/bin/env bash
# Import, build, and run all benchmark models

[[ -n "${BENCHMARK_COLLECTION}" ]] && model_dir="${BENCHMARK_COLLECTION}"

all_models=$(ls -1d ${BENCHMARK_DIR}/*/)

expected_to_work="
Boehm_JProteomeRes2014
Borghans_BiophysChem1997
Elowitz_Nature2000
Schwen_PONE2014
Chen_MSB2009
Fujita_SciSignal2010
Sneyd_PNAS2002
Zheng_PNAS2012"

for MODEL in $expected_to_work; do
    printf '=%.0s' {1..20}
    printf ${MODEL}
    printf '=%.0s' {1..20}
    echo
    echo ./import_and_run.sh ${BENCHMARK_DIR}/${MODEL}
    ./import_and_run.sh ${BENCHMARK_DIR}/${MODEL}
    printf '=%.0s' {1..100}
    echo
done
