#!/usr/bin/env bash
# Download and build Ipopt and HSL

set -euo pipefail
script_dir=$(dirname "$0")
script_dir=$(cd "${script_dir}" && pwd)
cd "${script_dir}"

make_opts=${MAKEOPTS-}

ipopt_url="https://www.coin-or.org/download/source/Ipopt/Ipopt-3.13.3.tgz"
ipopt_archive="Ipopt-3.13.3.tgz"
ipopt_dir="${script_dir}/Ipopt-releases-3.13.3"
ipopt_install_dir="${ipopt_dir}/install"
hsl_install_dir="${ipopt_dir}/ThirdParty-HSL/install"

if [[ ! -d "${ipopt_dir}" ]]; then
  if [[ ! -f "${ipopt_archive}" ]]; then
    echo "Downloading IpOpt source archive ..."
    wget -O "${ipopt_archive}" "${ipopt_url}"
  else
    echo "Skipping download step."
  fi

  tar -xzf "${ipopt_archive}"

  # Handle HSL solvers
  cd "${ipopt_dir}"
  if [[ ! -d "ThirdParty-HSL" ]]; then
    echo "Cloning ThirdParty-HSL ..."
    git clone https://github.com/coin-or-tools/ThirdParty-HSL.git
  fi

  cd ThirdParty-HSL

  # Need some coinhsl library.
  if [[ ! -d "coinhsl" ]]; then
    # ThirdParty/coinhsl ?
    if [[ -d "${script_dir}/coinhsl" ]]; then
      echo "Using coinhsl from ${script_dir}/coinhsl"
      cp -aR "${script_dir}/coinhsl" coinhsl
    else
      # Check for common ones:
      coinhsl_archive_names="
      coinhsl-2014.01.10
      coinhsl-2015.06.23
      coinhsl-2019.05.21
      "
      for coinhsl_archive in $coinhsl_archive_names; do
        if [[ -f "${script_dir}/${coinhsl_archive}.tar.gz" ]]; then
          echo "Unpacking ${script_dir}/${coinhsl_archive}.tar.gz ..."
          tar -xzf "${script_dir}/${coinhsl_archive}.tar.gz"
          mv "${coinhsl_archive}" coinhsl
        fi
      done
    fi

    if [[ ! -d "coinhsl" ]]; then
      echo "Did not find coinhsl/ or a known coinhsl archive." \
      "IpOpt will probably not work."
    fi
  fi

  # Use Intel MKL for lapack?
  lapack_lflags=""
  if [[ -v MKL_LIB ]]; then
    # will require F77=ifort when using intel compilers
    lapack_lflags="--with-lapack-lflags=${MKL_LIB}"
  fi

  ./configure --prefix="${hsl_install_dir}" --with-lapack "${lapack_lflags}" --enable-static --disable-shared
  make $make_opts
  make install

  cd "${ipopt_dir}"

  PKG_CONFIG_PATH=$PKG_CONFIG_PATH:${hsl_install_dir}/lib/pkgconfig/ \
  ./configure --prefix="${ipopt_install_dir}" \
    --enable-static \
    --disable-shared \
    --with-pic \
    --with-hsl \
    --disable-linear-solver-loader \
    --with-lapack "${lapack_lflags}"
    make $make_opts
  make $make_opts install
else
  echo "Skipping building Ipopt. Directory already exists."
fi
