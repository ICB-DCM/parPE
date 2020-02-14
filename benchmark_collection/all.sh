#!/usr/bin/env bash
# Import, build, and run benchmark models

set -e

[[ -n "${BENCHMARK_COLLECTION}" ]] && model_dir="${BENCHMARK_COLLECTION}"

# all_models=$(ls -1d ${model_dir}/*/)

expected_to_work="
Boehm_JProteomeRes2014
Borghans_BiophysChem1997
Elowitz_Nature2000
Chen_MSB2009
Fujita_SciSignal2010
Sneyd_PNAS2002
Zheng_PNAS2012"
# Schwen_PONE2014 TODO: shadowing observables


for MODEL in $expected_to_work; do
  printf '=%.0s' {1..20}
  printf %s "${MODEL}"
  printf '=%.0s' {1..20}
  echo
  echo ./import_and_run.sh "${model_dir}"/"${MODEL}"
  ./import_and_run.sh "${model_dir}"/"${MODEL}"
  printf '=%.0s' {1..100}
  echo
done
