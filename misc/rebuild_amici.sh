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
# delete old swig-generated files which can lead to ImportErrors during
# installation if the swig interface changed
set -x
rm -f "${amici_dir}/python/sdist/amici/amici_without_hdf5.py"
rm -f "${amici_dir}/python/sdist/amici/amici.py"
rm -f "${amici_dir}/python/sdist/amici/amici_wrap_without_hdf5.cxx"
rm -f "${amici_dir}/python/sdist/amici/amici_wrap.cxx"
rm -f "${amici_dir}"/python/sdist/amici/_amici.*.so
ls -l "${amici_dir}/python/sdist/amici"
# rebuild
"${amici_dir}/scripts/buildAmici.sh"

printf '=%.0s' {1..100}
echo "AMICI Python...."
# reinstall
"${script_dir}/run_in_venv.sh" "${build_dir}"/venv pip3 install -ve "${amici_dir}/python/sdist"


printf '=%.0s' {1..100}
