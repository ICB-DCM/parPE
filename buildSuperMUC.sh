#!/bin/bash
# Build parPE on SuperMUC

build_cpputest() {
  #cpputest
  #CPPUTEST_PATH=${amici_path}/ThirdParty/cpputest-master
  #cd ${CPPUTEST_PATH}
  # -DC++11=ON breaks compilation of some `_override` for no obvious reason
  #cmake -DC++11=OFF -DCMAKE_INSTALL_PREFIX=`pwd`
  #make ${make_opts}
  #make install
  :
}

build_amici() {
  amici_path=${parpe_root}/deps/AMICI

  "${amici_path}/scripts/buildSuiteSparse.sh"
  "${amici_path}/scripts/buildSundials.sh"

  python_exec=$(which python3)
  export PYTHON_EXECUTABLE=${python_exec}
  mkdir -p "${amici_path}/build"
  cd "${amici_path}/build"
  cmake -S .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DBLAS=MKL \
    -DBLAS_LIBRARIES="${MKL_LIB}" \
    -DBLAS_INCLUDE_DIRS="${MKL_INCDIR}" \
    -DBUILD_TESTS=OFF \
    -DBUILD_TESTING=OFF

  make ${make_opts} amici
}

build_3rd_party_deps() {
  # build dependencies
  cd "${parpe_root}/ThirdParty"
  #./installCeres.sh
  #./installCpputest.sh
  #  ceres_base=${parpe_root}/ThirdParty/ceres-solver-1.13.0/
  #  ceres_install_dir=${ceres_base}/build/install/
  #  if [[ -d ${ceres_base} ]]; then
  #    echo "Found CERES. Building..."
  #    mkdir -p "${ceres_base}/build"
  #    cd "${ceres_base}/build"
  #    make ${make_opts}
  #  else
  #    echo "CERES sources not found. Skipping..."
  #  fi
}

build_parpe() {
  echo
  echo "Building parPE..."
  cd "${parpe_root}"
  mkdir -p build && cd build
  HDF5_ROOT=${HDF5_BASE} \
  BOOST_ROOT=${BOOST_BASE} \
  MPI_HOME=${MPI_BASE} \
  cmake -S .. \
    -DBoost_USE_STATIC_LIBS=TRUE
  make ${make_opts}
}

set -eou pipefail

# absolute path to parpe repository base directory
parpe_root=$(dirname "$BASH_SOURCE")
parpe_root=$(cd "${parpe_root}" && pwd)

make_opts=${MAKEOPTS-}

build_cpputest
build_amici
build_3rd_party_deps
build_parpe
