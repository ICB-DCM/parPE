#!/bin/bash

export PARPE_BASE=`pwd`

# Build dependencies

# Install AMICI
AMICI_PATH=$PARPE_BASE/deps/AMICI/
cd $AMICI_PATH && scripts/buildSuiteSparse.sh && scripts/buildSundials.sh && scripts/buildCpputest.sh #&& scripts/buildAmici.sh
mkdir -p ${AMICI_PATH}/build && cd ${AMICI_PATH}/build
CPPUTEST_BUILD_DIR=${AMICI_PATH}/ThirdParty/cpputest-master/build/
cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_PYTHON=ON -DBUILD_TESTS=OFF -DCppUTest_DIR=${CPPUTEST_BUILD_DIR} .. && make -j12

#- cd $PARPE_BASE/ThirdParty && ./downloadPackages.sh
#- cd $PARPE_BASE/ThirdParty && ./installCeres.sh

cd $PARPE_BASE && ThirdParty/installGoogleTest.sh

# build parPE
pip install -r $PARPE_BASE/python/requirements.txt
cd $PARPE_BASE && ./buildShippable.sh

# run tests
cd $PARPE_BASE/build && CTEST_OUTPUT_ON_FAILURE=1 make test

# valgrind
CTEST_OUTPUT_ON_FAILURE=1 make ExperimentalMemCheck; cat Testing/Temporary/MemoryChecker.*.log
