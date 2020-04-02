#!/bin/bash
# build parPE dependencies

script_dir=$(dirname "$0")
script_dir=$(cd "$script_dir" && pwd )

set -euo pipefail

cd "$script_dir"

../deps/AMICI/scripts/buildAll.sh

./installCBLAS.sh

./installCeres.sh

./installIpopt.sh
