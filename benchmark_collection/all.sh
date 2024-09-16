#!/usr/bin/env bash
# Import, build, and run benchmark models

set -eou pipefail

# Path to benchmark collection model directory
[[ -n "${BENCHMARK_COLLECTION}" ]] && model_dir="${BENCHMARK_COLLECTION}"

# all_models=$(ls -1d ${model_dir}/*/)
expected_to_work="
Boehm_JProteomeRes2014
Borghans_BiophysChem1997
Crauste_CellSystems2017
Elowitz_Nature2000
Fujita_SciSignal2010
Sneyd_PNAS2002
Weber_BMC2015
Zheng_PNAS2012
"
# Timepoint-specific parameter overrides currently unsupported.:
#
# Chen_MSB2009
#
# missing ref:
#
# Alkan_SciSignal2018
# Blasi_CellSystems2016
# Bruno_JExpBio2016
# Giordano_Nature2020
# Okuonghae_ChaosSolitonsFractals2020
# Perelson_Science1996
# Rahman_MBS2016
# Raimundez_PCB2020
# SalazarCavazos_MBoC2020
# Zhao_QuantBiol2020
#
# integration trouble:
#
# Beer_MolBioSystems2014
# Bertozzi_PNAS2020
# Borghans_BiophysChem1997
# Isensee_JCB2018
#
# Mismatch:
# Bachmann_MSB2011: Expected llh -478.459689232875, got nllh -418.409381
# Brannmark_JBC2010 FAILED: Expected llh 283.778227541074, got nllh 141.889034
# Lucarelli_CellSystems2018
# Fiedler_BMC2016 FAILED: Expected llh -117.16780323362, got nllh 109391.496205
#
# Unsupported priors:
#   Schwen_PONE2014
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
