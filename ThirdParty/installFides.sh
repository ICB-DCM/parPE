#!/usr/bin/env bash
# Download and build fides-cpp

set -euo pipefail
script_dir=$(dirname "$0")
script_dir=$(cd "${script_dir}" && pwd)
cd "${script_dir}"

make_opts=${MAKEOPTS-}

fides_dir=${script_dir}/fides-cpp

if [ ! -d "${fides_dir}" ]; then
  cd "${script_dir}" && git clone https://github.com/dweindl/fides-cpp.git
fi

cd "${fides_dir}"
cmake -S . -B build
cd build
make ${make_opts} fides
