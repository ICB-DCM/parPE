#!/usr/bin/env bash
# Install Ceres solver and Eigen lib

set -euo pipefail
set -x
script_dir=$(dirname "$0")
script_dir=$(cd "${script_dir}" && pwd )

make_opts=${MAKEOPTS-}
cd "${script_dir}"

# build EIGEN
eigen_version="3.3.9"
eigen_url="https://gitlab.com/libeigen/eigen/-/archive/${eigen_version}/eigen-${eigen_version}.tar.gz"
eigen_archive="eigen-${eigen_version}.tar.gz"
eigen_dir="${script_dir}/eigen-${eigen_version}"
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
ceres_version="2.0.0"
ceres_url="http://ceres-solver.org/ceres-solver-${ceres_version}.tar.gz"
ceres_archive="ceres-solver-${ceres_version}.tar.gz"
ceres_dir="${script_dir}/ceres-solver-${ceres_version}"
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
        -DMINIGLOG=ON \
        ..

  make $make_opts
  make $make_opts install
fi
