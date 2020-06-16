#!/usr/bin/env bash
# Rebuild AMICI and dependencies

set -euo pipefail

script_dir=$(dirname "$0")
script_dir=$(cd "$script_dir" && pwd)

build_dir="${script_dir}/../build"
amici_dir="${script_dir}/../deps/AMICI"

printf '=%.0s' {1..100}
echo "Suitesparse...."
"${amici_dir}/scripts/buildSuiteSparse.sh"

printf '=%.0s' {1..100}
echo "Sundials...."
"${amici_dir}/scripts/buildSundials.sh"

printf '=%.0s' {1..100}
echo "AMICI C++...."
"${amici_dir}/scripts/buildAmici.sh"

printf '=%.0s' {1..100}
echo "AMICI Python...."
"${script_dir}/run_in_venv.sh" "${build_dir}"/venv pip3 install -ve "${amici_dir}/python/sdist"


printf '=%.0s' {1..100}
