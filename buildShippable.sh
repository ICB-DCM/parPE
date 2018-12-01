#!/bin/bash
# Build script for CI on shippable.com

mkdir -p build
cd build

CC=mpicc CXX=mpiCC cmake \
      -DIPOPT_INCLUDE_DIRS=/usr/include/coin/ \
      -DIPOPT_LIBRARIES=/usr/lib/libipopt.so \
      -DCERES_LIBRARIES="/usr/lib/libceres.so;/usr/lib/x86_64-linux-gnu/libglog.so;/usr/lib/x86_64-linux-gnu/libgflags.so" \
      -DCERES_INCLUDE_DIRS="/usr/include/;/usr/include/eigen3" \
      -DCppUTest_ROOT=`pwd`/../deps/AMICI/ThirdParty/cpputest-master/build/ \
      -DMPI_INCLUDE_DIRS=/usr/include/openmpi-x86_64/ \
      -DGCOVR_REPORT=TRUE \
      ..

make -j12 VERBOSE=1

cd ..
