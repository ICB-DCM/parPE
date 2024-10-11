#!/bin/bash
# Install parPE inside docker
set -euo pipefail
set -x

export PARPE_BASE="${PARPE_DIR:-/parPE}"

# unpack git archive
mkdir "$PARPE_BASE" && cd "$PARPE_BASE"
tar -xzf /container-files/parpe.tar.gz

# Build dependencies

# Install AMICI
export AMICI_PATH=${PARPE_BASE}/deps/AMICI/
cd "${AMICI_PATH}" \
  && scripts/buildSuiteSparse.sh \
  && scripts/buildSundials.sh
mkdir -p "${AMICI_PATH}/build" && cd "${AMICI_PATH}/build"
cmake \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DENABLE_PYTHON=ON \
  -DBUILD_TESTS=OFF \
  .. && make -j

# install fides optimizer
cd "$PARPE_BASE/ThirdParty" && ./installFides.sh

# build parPE
cd "${PARPE_BASE}"
mkdir -p build
cd build
ceres_libs="/usr/lib/libceres.so"
ceres_libs="$ceres_libs;/usr/lib/x86_64-linux-gnu/libglog.so"
ceres_libs="$ceres_libs;/usr/lib/x86_64-linux-gnu/libgflags.so"

# docker-specific MPI settings
mpi_cmd="mpiexec;--allow-run-as-root;-n;4;--oversubscribe"
mpi_cmd="$mpi_cmd;--mca;btl_vader_single_copy_mechanism;none"
mpi_cmd="$mpi_cmd;--mca;btl;^openib"
mpi_cmd="$mpi_cmd;--mca;oob_tcp_if_include;lo"
mpi_cmd="$mpi_cmd;--mca;btl_tcp_if_include;lo;"
mpi_cmd="$mpi_cmd;--mca;orte_base_help_aggregate;0"

CC=mpicc CXX=mpiCC cmake \
      -DPARPE_BUILD_OPTIMIZED=OFF \
      -DPARPE_ENABLE_FIDES=ON \
      -DIPOPT_INCLUDE_DIRS=/usr/include/coin/ \
      -DIPOPT_LIBRARIES=/usr/lib/libipopt.so \
      -DMPI_INCLUDE_DIRS=/usr/include/openmpi-x86_64/ \
      -DBUILD_TESTING=ON \
      -DTESTS_MPIEXEC_COMMAND="$mpi_cmd" \
      ..

make -j VERBOSE=1

# MPI settings for python tests
export PARPE_TESTS_MPIEXEC="mpiexec -n 5 --oversubscribe --allow-run-as-root --mca btl_vader_single_copy_mechanism none --mca btl ^openib --mca oob_tcp_if_include lo --mca btl_tcp_if_include lo --mca orte_base_help_aggregate 0"

# run tests
cd "${PARPE_BASE}"/build && CTEST_OUTPUT_ON_FAILURE=1 make test

# valgrind
#CTEST_OUTPUT_ON_FAILURE=1 make ExperimentalMemCheck; cat Testing/Temporary/MemoryChecker.*.log
