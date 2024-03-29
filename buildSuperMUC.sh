#!/bin/bash
# Build parPE on SuperMUC
#
# Assumes all relevant modules are loaded with compatible versions.
# Assumes there is Python>=3.7 available on $PATH.
#
# Known issues:
# * With CMake<3.16, linking IpOpt compiled with MKL via pkg-config info may
#   result in "ld: cannot find -lpkgcfg_lib_IPOPT_iomp5-NOTFOUND"
#   -> use newer CMake

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
    -DBUILD_TESTS=OFF
  make ${make_opts} amici
}

build_boost() {
  # Build custom boost (with CXX11 ABI!)
  cd "${parpe_root}/ThirdParty"
  # NOTE 1.70.0 has broken CMake config
  local archive_name=boost_1_71_0.tar.gz
  local boost_dir=${parpe_root}/ThirdParty/boost_1_71_0
  local url="https://dl.bintray.com/boostorg/release/1.71.0/source/${archive_name}"
  boost_install_dir=${boost_dir}/install

  if [[ ! -d ${boost_dir} ]]; then
    if [[ ! -f ${archive_name} ]];
      then wget $url -O "${archive_name}"
    fi

    tar -xzf $archive_name
    cd "${boost_dir}"
    ./bootstrap.sh \
      --prefix="${boost_install_dir}" \
      --with-libraries=filesystem,program_options,regex,serialization,system
    ./b2 install -j 20
  else
    echo "Skipping boost - ${boost_dir} already exists"
  fi
}

build_3rd_party_deps() {
  # build dependencies
  # build_boost
  "${parpe_root}/ThirdParty/installIpopt.sh"
  #./installCeres.sh
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
  ipopt_root=${parpe_root}/ThirdParty/Ipopt-releases-3.13.3/

  #BOOST_ROOT=${boost_install_dir} \
  HDF5_ROOT=${HDF5_BASE} \
  MPI_HOME=${MPI_BASE} \
  PKG_CONFIG_PATH=${PKG_CONFIG_PATH:-}:${ipopt_root}/install/lib/pkgconfig/:${ipopt_root}/ThirdParty-HSL/install/lib/pkgconfig/ \
  cmake -S .. \
    -DPARPE_ENABLE_CERES=OFF \
    -DBoost_USE_STATIC_LIBS=TRUE \
    -DBUILD_TESTING=OFF
  make ${make_opts}
}

set -eou pipefail

# absolute path to parpe repository base directory
parpe_root=$(dirname "$BASH_SOURCE")
parpe_root=$(cd "${parpe_root}" && pwd)

make_opts=${MAKEOPTS--j12}

build_3rd_party_deps
build_amici
build_parpe
