#!/bin/bash

set -euo pipefail

parpe_root=$(dirname "$0")
parpe_root=$(cd "${parpe_root}" && pwd)

# Build amici
amici_dir="${parpe_root}/deps/AMICI/"
cd "${amici_dir}"
#scripts/buildAll.sh

# build ceres
"${parpe_root}/ThirdParty/installCeres.sh"
# build ipopt
"${parpe_root}/ThirdParty/installIpopt.sh"

# build parpe
ceres_install_dir="${parpe_root}/ThirdParty/ceres-solver-1.13.0/build/install"
ceres_libs="${ceres_install_dir}/lib/libceres.a;$(pkg-config --libs blas);cxsparse;lapack;cholmod;camd;colamd"
ceres_inc="${ceres_install_dir}/include/"
ceres_inc="${ceres_inc};${parpe_root}/ThirdParty/eigen-eigen-67e894c6cd8f/build/install/include/eigen3/"
parpe_build_dir="${parpe_root}/build"
mkdir -p "${parpe_build_dir}"
cd "${parpe_build_dir}"
CC=mpicc CXX=mpiCC cmake \
      -DIPOPT_DIR="${parpe_root}/ThirdParty/Ipopt-3.12.12/install" \
      -DCERES_LIBRARIES="${ceres_libs}" \
      -DCERES_INCLUDE_DIRS="${ceres_inc}" \
      -DENABLE_SWIG=FALSE \
      "${parpe_root}"
make -j12

# run tests
make test
