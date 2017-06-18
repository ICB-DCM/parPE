#!/bin/bash
mkdir build
cd build
CC=mpicc CXX=mpiCC cmake -DAMICI_DIR=`pwd`/../ThirdParty/AMICI-master \
      -DIPOPT_INCLUDE_DIRS=/usr/include/coin/ \
      -DIPOPT_LIBRARIES=/usr/lib/libipopt.a \
      -DCERES_LIBRARIES="/usr/lib/libceres.a;/usr/lib/x86_64-linux-gnu/libglog.a;/usr/lib/x86_64-linux-gnu/libgflags.a" \
      -DCERES_INCLUDE_DIRS="/usr/include/;/usr/include/eigen3" \
      -DCPPUTEST_DIR=`pwd`/../ThirdParty/AMICI-master/ThirdParty/cpputest-3.8/ \
      -DMPI_INCLUDE_DIR=/usr/include/openmpi-x86_64/ \
      -DMPI_LIBRARY=/usr/lib64/openmpi/lib/libmpi_cxx.so \
      -DBLAS_INCLUDE_DIRS="/usr/include/" \
      -DBLAS_LIBRARIES="/usr/lib/libcblas.a;/usr/lib/libatlas.a;/usr/lib/libblas.a" \
      ..
make -j12 VERBOSE=1
