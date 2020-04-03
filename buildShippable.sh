#!/bin/bash
# Build script for CI on shippable.com
set -euo pipefail
set -x

# build parPE
build_dir="$(pwd)/build"
CC=mpicc CXX=mpiCC cmake \
      -B"${build_dir}" \
      -DIPOPT_INCLUDE_DIRS=/usr/include/coin/ \
      -DIPOPT_LIBRARIES=/usr/lib/libipopt.so \
      -DCERES_LIBRARIES="/usr/lib/libceres.so;/usr/lib/x86_64-linux-gnu/libglog.so;/usr/lib/x86_64-linux-gnu/libgflags.so" \
      -DCERES_INCLUDE_DIRS="/usr/include/;/usr/include/eigen3" \
      -DMPI_INCLUDE_DIRS=/usr/include/openmpi-x86_64/ \
      -DGCOVR_REPORT=TRUE \
      -DBUILD_TESTS=TRUE \
      ..

build-wrapper-linux-x86-64 --out-dir bw-output cmake --build "${build_dir}" -j12 -- VERBOSE=1

sonar-scanner \
  -Dsonar.organization=icb-dcm \
  -Dsonar.projectKey=ICB-DCM_parPE \
  -Dsonar.sources=. \
  -Dsonar.host.url=https://sonarcloud.io \
  -Dsonar.cfamily.build-wrapper-output=bw-output \
  -Dsonar.login=45879b85332d963bf3ead20399bf5bd2c925c156
