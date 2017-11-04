#!/bin/bash
set -e 

SCRIPT_PATH="`dirname \"$0\"`"
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"

cd ${SCRIPT_PATH}/../deps/AMICI/ThirdParty/cpputest-master
mkdir -p build-noleakcheck
cd build-noleakcheck
cmake -DC++11=ON -DBUILD_TESTING=OFF -DMEMORY_LEAK_DETECTION=OFF -DCMAKE_INSTALL_PREFIX=`pwd`/install ..
make -j12 install
