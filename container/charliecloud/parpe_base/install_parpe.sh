#!/bin/bash
set -e 
cd
git clone --single-branch --branch feature_charlie --depth=1 https://github.com/ICB-DCM/parPE.git
cd parPE
export PARPE_BASE=`pwd`

# FIXME
sed -ri 's/"README\.md", "r"/"README\.md", "r", encoding="utf-8"/' deps/AMICI/python/sdist/setup.py

# Build dependencies

# Install AMICI
export AMICI_PATH=${PARPE_BASE}/deps/AMICI/
cd $AMICI_PATH && scripts/buildSuiteSparse.sh && scripts/buildSundials.sh && scripts/buildCpputest.sh #&& scripts/buildAmici.sh
mkdir -p ${AMICI_PATH}/build && cd ${AMICI_PATH}/build
CPPUTEST_BUILD_DIR=${AMICI_PATH}/ThirdParty/cpputest-master/build/
cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_PYTHON=ON -DBUILD_TESTS=OFF -DCppUTest_DIR=${CPPUTEST_BUILD_DIR} .. && make -j12

#- cd $PARPE_BASE/ThirdParty && ./downloadPackages.sh
#- cd $PARPE_BASE/ThirdParty && ./installCeres.sh

cd $PARPE_BASE && ThirdParty/installGoogleTest.sh

# build parPE
pip install -r $PARPE_BASE/python/requirements.txt
cd $PARPE_BASE
mkdir -p build
cd build
CC=mpicc CXX=mpiCC cmake \
      -DIPOPT_INCLUDE_DIRS=/usr/include/coin/ \
      -DIPOPT_LIBRARIES=/usr/lib/libipopt.so \
      -DCERES_LIBRARIES="/usr/lib/libceres.so;/usr/lib/x86_64-linux-gnu/libglog.so;/usr/lib/x86_64-linux-gnu/libgflags.so" \
      -DCERES_INCLUDE_DIRS="/usr/include/;/usr/include/eigen3" \
      -DMPI_INCLUDE_DIRS=/usr/include/openmpi-x86_64/ \
      ..
make -j12 VERBOSE=1


# run tests
cd $PARPE_BASE/build && CTEST_OUTPUT_ON_FAILURE=1 make test

# valgrind
CTEST_OUTPUT_ON_FAILURE=1 make ExperimentalMemCheck; cat Testing/Temporary/MemoryChecker.*.log
