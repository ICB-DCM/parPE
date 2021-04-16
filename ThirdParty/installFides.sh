#!/usr/bin/env bash
# Download and build fides-cpp

set -euo pipefail
script_dir=$(dirname "$0")
script_dir=$(cd "${script_dir}" && pwd)
cd "${script_dir}"

make_opts=${MAKEOPTS-}

fides_dir=${script_dir}/fides-cpp
blaze_dir=${script_dir}/blaze-3.8

cd "${script_dir}"

if [[ ! -d "${blaze_dir}" ]]; then
  if [[ ! -f "blaze-3.8.tar.gz" ]]; then
    wget https://bitbucket.org/blaze-lib/blaze/downloads/blaze-3.8.tar.gz
  fi
  tar -xzf blaze-3.8.tar.gz
  cd blaze-3.8
  cmake -S . -B build/ -DCMAKE_INSTALL_PREFIX=/usr/local/
  cd build
  sudo make -j2 install
fi

if [ ! -d "${fides_dir}" ]; then
  cd "${script_dir}" && git clone https://github.com/dweindl/fides-cpp.git
fi

cd "${fides_dir}"
cmake -S . -B build
cd build
make ${make_opts}
