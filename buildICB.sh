#!/bin/bash
mkdir build
cd build
CC=mpicc CXX=mpiCC cmake -DAMICI_DIR=$HOME/src/AMICI \
      -DIPOPT_DIR=`pwd`/../ThirdParty/Ipopt-3.12.7/install \
      -DCERES_LIBRARIES=`pwd`/../ThirdParty/ceres-solver-1.13.0/build/install/lib64/libceres.a \
      -DCERES_INCLUDE_DIRS="`pwd`/../ThirdParty/ceres-solver-1.12.0/build/install/include/;`pwd`/../ThirdParty/ceres-solver-1.12.0/build/install/include/ceres/internal/miniglog/;`pwd`/../ThirdParty/eigen-eigen-67e894c6cd8f/build/install/include/eigen3/" \
      -DCPPUTEST_DIR=`pwd`/../ThirdParty/cpputest-3.8/ \
      -DMPI_INCLUDE_DIR=/usr/include/openmpi-x86_64/ \
      -DMPI_LIBRARY=/usr/lib64/openmpi/lib/libmpi_cxx.so \
      -DBLAS_INCLUDE_DIRS=`pwd`/../ThirdParty/CBLAS/include \
      -DBLAS_LIBRARIES=`pwd`/../ThirdParty/CBLAS/lib/cblas_LINUX.a \
      ..
make -j12
