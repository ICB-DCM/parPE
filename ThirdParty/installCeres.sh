#!/usr/bin/env bash
# Install Ceres solver and Eigen lib

set -euo pipefail
set -x
script_dir=$(dirname "$0")
script_dir=$(cd "${script_dir}" && pwd )

make_opts=${MAKEOPTS-}
cd "${script_dir}"

# build EIGEN
eigen_url="http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz"
eigen_archive="eigen-3.3.3.tar.gz"
eigen_dir="${script_dir}/eigen-eigen-67e894c6cd8f"
eigen_build_dir="${eigen_dir}/build"
eigen_install_dir="${eigen_build_dir}/install"

if [[ ! -d "${eigen_dir}" ]]; then
  if [[ ! -f "${eigen_archive}" ]]
      then wget -O "${eigen_archive}" "${eigen_url}"
  fi

  tar -xzf "${eigen_archive}"

  mkdir -p "${eigen_build_dir}"
  cd "${eigen_build_dir}"
  cmake -DCMAKE_INSTALL_PREFIX="$eigen_install_dir" ..
  make ${make_opts}
  make ${make_opts} install
fi


# build CERES
cd "$script_dir"
ceres_url="http://ceres-solver.org/ceres-solver-1.13.0.tar.gz"
ceres_archive="ceres-solver-1.13.0.tar.gz"
ceres_dir="${script_dir}/ceres-solver-1.13.0"
ceres_build_dir="${ceres_dir}/build"
ceres_install_dir="${ceres_build_dir}/install"

if [[ ! -d "${ceres_dir}" ]]; then
  if [[ ! -f "${ceres_archive}" ]]
      then wget -O "${ceres_archive}" "${ceres_url}"
  fi

  tar -xzf "${ceres_archive}"

  mkdir -p "${ceres_build_dir}"
  cd "${ceres_build_dir}"
  cmake -DCMAKE_INSTALL_PREFIX="${ceres_install_dir}" \
        -DBUILD_SHARED_LIBS=OFF \
        -DBUILD_TESTING=OFF \
        -DBUILD_EXAMPLES=OFF \
        -DLAPACK=OFF \
        -DGFLAGS=OFF \
        -DCXX11=ON \
        -DEIGEN_INCLUDE_DIR="${eigen_dir}" \
        -DEIGEN_INCLUDE_DIR_HINTS="${eigen_install_dir}/include" \
        -DMINIGLOG=ON \
        ..

  make $make_opts
  make $make_opts install
fi
