#!/bin/bash
set -e
PARPE_ROOT="`dirname \"$BASH_SOURCE\"`"
PARPE_ROOT="`( cd \"${PARPE_ROOT}\" && pwd )`"

MAKE_OPTS=-j12
# load modules

# build AMICI
AMICI_PATH=${PARPE_ROOT}/deps/AMICI

#cpputest
CPPUTEST_PATH=${AMICI_PATH}/ThirdParty/cpputest-master
cd ${CPPUTEST_PATH}
# -DC++11=ON breaks compilation of some `_override` for no obvious reason
cmake -DC++11=OFF -DCMAKE_INSTALL_PREFIX=`pwd` 
make ${MAKE_OPTS}
make install

mkdir -p ${AMICI_PATH}/build
cd ${AMICI_PATH}/build
cmake -DCMAKE_BUILD_TYPE=Debug \
	-DBoost_INCLUDE_DIR=${BOOST_INCDIR} \
	-DBLAS=MKL \
	-DBLAS_LIBRARIES=${MKL_LIB} \
	-DBLAS_INCLUDE_DIRS=${MKL_INCDIR} \
	 ..
make ${MAKE_OPTS}

# build dependencies
cd ${PARPE_ROOT}/ThirdParty
#./installCeres.sh
#./installCpputest.sh
# requires download of additional packages ./installIpopt.sh

CERES_BASE=${PARPE_ROOT}/ThirdParty/ceres-solver-1.13.0/
mkdir -p ${CERES_BASE}/build
cd ${CERES_BASE}/build
CC=mpicc CXX=mpiCC HDF5_ROOT=${HDF5_BASE} BOOSTROOT=${BOOST_BASE} MPI_HOME=${MPI_BASE} cmake \
      -DIPOPT_DIR=`pwd`/../ThirdParty/Ipopt-3.12.7/install \
      -DCERES_LIBRARIES=${CERES_BASE}/build/install/lib64/libceres.a \
      -DCERES_INCLUDE_DIRS="${CERES_BASE}/build/install/include/;${CERES_BASE}/build/install/include/ceres/internal/miniglog/;`pwd`/../ThirdParty/eigen-eigen-67e894c6cd8f/build/install/include/eigen3/" \
      -DCPPUTEST_DIR=${CPPUTEST_PATH} \
      ..
make ${MAKE_OPTS}
