#!/bin/bash

set -e 

# build EIGEN
tar -xzf eigen-3.3.3.tar.gz
cd eigen-eigen-67e894c6cd8f
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=`pwd`/install ..
make -j12
make install
cd ../..

# build CERES
tar -xzf ceres-solver-1.13.0.tar.gz
cd ceres-solver-1.13.0/
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=`pwd`/install \
      -DBUILD_SHARED_LIBS=ON \
      -DBUILD_TESTING=OFF \
      -DBUILD_EXAMPLES=OFF \
      -DGFLAGS=OFF \
      -DCXX11=ON \
      -DEIGEN_INCLUDE_DIR=`pwd`/../../eigen-eigen-67e894c6cd8f/ \
      -DMINIGLOG=ON \
      ..
      
make -j12
make install

      
