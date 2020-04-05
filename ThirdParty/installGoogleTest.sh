#!/usr/bin/env bash
# Download and build googletest

set -euo pipefail

script_path="$(dirname "$0")"
script_path="$(cd "${script_path}" && pwd)"

cd "${script_path}"

if [[ ! -d "googletest" ]]; then
  echo "googletest/ does not exist"
  if [[ ! -f "googletest.zip" ]]; then
    echo "Downloading googletest..."
    wget "https://github.com/google/googletest/archive/master.zip" -O "googletest.zip"
  fi

  echo "Unpacking and building googletest..."
  unzip "googletest"
  mv "googletest-master" "googletest"
  cd "googletest"
  mkdir -p build
  cd build
  cmake ..
  make -j 4
else
  echo "googletest/ exists. nothing to do."
fi
