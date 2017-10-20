#!/bin/bash

# load modules

# build dependencies
cd ThirdParty
#./installCeres.sh
#./installCpputest.sh
# requires download of additional packages ./installIpopt.sh
cd ..

CERES_BASE=`pwd`/ThirdParty/ceres-solver-1.13.0/

mkdir -p build
cd build
CC=mpicc CXX=mpiCC HDF5_ROOT=${HDF5_BASE} BOOSTROOT=${BOOST_BASE} MPI_HOME=${MPI_BASE} cmake \
      -DAMICI_DIR=$HOME/src/AMICI \
      -DBLAS=MKL \
      -DIPOPT_DIR=$HOME/src/Ipopt-3.12.6/install \
      -DCERES_LIBRARIES=${CERES_BASE}/build/install/lib64/libceres.so \
      -DCERES_INCLUDE_DIRS="${CERES_BASE}/build/install/include/;${CERES_BASE}/build/install/include/ceres/internal/miniglog/;`pwd`/../ThirdParty/eigen-eigen-67e894c6cd8f/build/install/include/eigen3/" \
      -DCPPUTEST_DIR=`pwd`/../ThirdParty/cpputest-3.8/ \
      -DBLAS_INCLUDE_DIRS="${MKL_INCDIR}" \
      -DBLAS_LIBRARIES="${MKL_LIB}" \
      ..
make -j12
cd ..
