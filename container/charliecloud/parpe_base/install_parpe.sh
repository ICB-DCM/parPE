#!/bin/bash
set -e

cd

# unpack git archive
mkdir parPE && cd parPE
tar -xzf /u18/parpe.tar.gz

export PARPE_BASE=$(pwd)

# Build dependencies

# Install AMICI
export AMICI_PATH=${PARPE_BASE}/deps/AMICI/
cd "${AMICI_PATH}" && scripts/buildSuiteSparse.sh && scripts/buildSundials.sh && scripts/buildCpputest.sh #&& scripts/buildAmici.sh
mkdir -p "${AMICI_PATH}"/build && cd "${AMICI_PATH}"/build
CPPUTEST_BUILD_DIR=${AMICI_PATH}/ThirdParty/cpputest-master/build/
cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_PYTHON=ON -DBUILD_TESTS=OFF -DCppUTest_DIR="${CPPUTEST_BUILD_DIR}" .. && make -j12

#- cd $PARPE_BASE/ThirdParty && ./downloadPackages.sh
#- cd $PARPE_BASE/ThirdParty && ./installCeres.sh

# For google-test for parPE tests
cd "${PARPE_BASE}" && ThirdParty/installGoogleTest.sh

# build parPE
pip install -r "${PARPE_BASE}"/python/requirements.txt
cd "${PARPE_BASE}"
mkdir -p build
cd build
CC=mpicc CXX=mpiCC cmake \
      -DIPOPT_INCLUDE_DIRS=/usr/include/coin/ \
      -DIPOPT_LIBRARIES=/usr/lib/libipopt.so \
      -DCERES_LIBRARIES="/usr/lib/libceres.so;/usr/lib/x86_64-linux-gnu/libglog.so;/usr/lib/x86_64-linux-gnu/libgflags.so" \
      -DCERES_INCLUDE_DIRS="/usr/include/;/usr/include/eigen3" \
      -DMPI_INCLUDE_DIRS=/usr/include/openmpi-x86_64/ \
      -DBUILD_TESTS=ON \
      "-DTESTS_MPIEXEC_COMMAND=mpiexec;--allow-run-as-root;-n;4;--oversubscribe;--mca;btl_vader_single_copy_mechanism;none" \
      ..
make -j12 VERBOSE=1

# run tests
cd "${PARPE_BASE}"/build && CTEST_OUTPUT_ON_FAILURE=1 make test

# valgrind
#CTEST_OUTPUT_ON_FAILURE=1 make ExperimentalMemCheck; cat Testing/Temporary/MemoryChecker.*.log
