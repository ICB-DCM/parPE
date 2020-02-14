#!/bin/bash
# build parPE dependencies

script_dir=$(dirname "$0")
script_dir=$(cd "$script_dir" && pwd )

# exit on error
set -e

../deps/AMICI/scripts/buildAll.sh

./installCBLAS.sh

./installCeres.sh

./installIpopt.sh
