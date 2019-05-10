#!/bin/bash
set -e
PARPE_ROOT="`dirname \"$BASH_SOURCE\"`"
PARPE_ROOT="`( cd \"${PARPE_ROOT}\" && pwd )`"

MAKE_OPTS=-j12
# load modules

# build AMICI
AMICI_PATH=${PARPE_ROOT}/deps/AMICI

#cpputest
#CPPUTEST_PATH=${AMICI_PATH}/ThirdParty/cpputest-master
#cd ${CPPUTEST_PATH}
# -DC++11=ON breaks compilation of some `_override` for no obvious reason
#cmake -DC++11=OFF -DCMAKE_INSTALL_PREFIX=`pwd`
#make ${MAKE_OPTS}
#make install

deps/AMICI/scripts/buildSuiteSparse.sh
deps/AMICI/scripts/buildSundials.sh

export PYTHON_EXECUTABLE=$(which python3)
mkdir -p ${AMICI_PATH}/build
cd ${AMICI_PATH}/build
cmake -DCMAKE_BUILD_TYPE=Release \
        -DBLAS=MKL \
        -DBLAS_LIBRARIES="${MKL_LIB}" \
        -DBLAS_INCLUDE_DIRS="${MKL_INCDIR}" \
	-DBUILD_TESTS=OFF \
        ..
make ${MAKE_OPTS} amici

# build dependencies
cd ${PARPE_ROOT}/ThirdParty
#./installCeres.sh
#./installCpputest.sh
# requires download of additional packages ./installIpopt.sh

CERES_BASE=${PARPE_ROOT}/ThirdParty/ceres-solver-1.13.0/
CERES_INSTALL_DIR=${CERES_BASE}/build/install/
if [[ -d ${CERES_BASE} ]]; then
    echo "Found CERES. Building..."
    mkdir -p ${CERES_BASE}/build
    cd ${CERES_BASE}/build
    make ${MAKE_OPTS}
else
   echo "CERES sources not found. Skipping..."
fi

echo
echo "Building parPE..."
cd $PARPE_ROOT
mkdir -p build && cd build
CC=mpicc CXX=mpiCC HDF5_ROOT=${HDF5_BASE} BOOST_ROOT=${BOOST_BASE} MPI_HOME=${MPI_BASE} cmake \
      -DBoost_USE_STATIC_LIBS=TRUE \
      -DIPOPT_DIR=`pwd`/../ThirdParty/Ipopt-3.12.9/install \
      -DCERES_LIBRARIES="${CERES_INSTALL_DIR}/lib64/libceres.a;${MKL_LIB}" \
      -DCERES_INCLUDE_DIRS="${CERES_INSTALL_DIR}/include/;${CERES_INSTALL_DIR}/include/ceres/internal/miniglog/;`pwd`/../ThirdParty/eigen-eigen-67e894c6cd8f/build/install/include/eigen3/" \
      ..
make ${MAKE_OPTS}
