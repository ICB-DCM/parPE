#!/usr/bin/env bash
# Import, build, and run benchmark models

set -eou pipefail

# Path to benchmark collection model directory
[[ -n "${BENCHMARK_COLLECTION}" ]] && model_dir="${BENCHMARK_COLLECTION}"

# all_models=$(ls -1d ${model_dir}/*/)

expected_to_work="
Boehm_JProteomeRes2014
Borghans_BiophysChem1997
Elowitz_Nature2000
Sneyd_PNAS2002
Zheng_PNAS2012"
# Fujita_SciSignal2010 "Timepoint-specific parameter overrides currently unsupported." in PEtab parameter mapping
# Schwen_PONE2014 Chen_MSB2009

for model_name in $expected_to_work; do
  printf '=%.0s' {1..20}
  printf %s "${model_name}"
  printf '=%.0s' {1..20}
  echo
  echo ./import_and_run.sh "${model_dir}"/"${model_name}"
  ./import_and_run.sh "${model_dir}"/"${model_name}"
  printf '=%.0s' {1..100}
  echo
done
