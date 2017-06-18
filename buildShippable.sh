#!/bin/bash
mkdir build
cd build
CC=mpicc CXX=mpiCC cmake -DAMICI_DIR=`pwd`/../ThirdParty/AMICI-master \
      -DIPOPT_INCLUDE_DIRS=/usr/include/coin/ \
      -DIPOPT_LIBRARIES=/usr/lib/libipopt.a \
      -DCERES_LIBRARIES=/usr/lib/libceres.a \
      -DCERES_INCLUDE_DIRS="/usr/include/" \
      -DCPPUTEST_DIR=`pwd`/../ThirdParty/AMICI-master/ThirdParty/cpputest-3.8/ \
      -DMPI_INCLUDE_DIR=/usr/include/openmpi-x86_64/ \
      -DMPI_LIBRARY=/usr/lib64/openmpi/lib/libmpi_cxx.so \
      -DBLAS_INCLUDE_DIRS="/usr/include/" \
      -DBLAS_LIBRARIES=/usr/lib/libcblas.a \
      ..
make -j12 VERBOSE=1
