#!/bin/bash

mkdir build
cd build

CC=mpicc CXX=mpiCC cmake \
      -DAMICI_DIR=`pwd`/../ThirdParty/AMICI-master \
      -DIPOPT_INCLUDE_DIRS=/usr/include/coin/ \
      -DIPOPT_LIBRARIES=/usr/lib/libipopt.so \
      -DCERES_LIBRARIES="/usr/lib/libceres.so;/usr/lib/x86_64-linux-gnu/libglog.so;/usr/lib/x86_64-linux-gnu/libgflags.so" \
      -DCERES_INCLUDE_DIRS="/usr/include/;/usr/include/eigen3" \
      -DCPPUTEST_DIR=`pwd`/../ThirdParty/AMICI-master/ThirdParty/cpputest-3.8/ \
      -DMPI_INCLUDE_DIR=/usr/include/openmpi-x86_64/ \
      -DMPI_LIBRARY=/usr/lib64/openmpi/lib/libmpi_cxx.so \
      -DBLAS_INCLUDE_DIRS="/usr/include/" \
      -DBLAS_LIBRARIES="/usr/lib/libcblas.so;/usr/lib/libatlas.so;/usr/lib/libblas.so" \
      -DGCOVR_REPORT=TRUE \
      ..

make -j12 VERBOSE=1

cd ..
