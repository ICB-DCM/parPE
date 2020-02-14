#!/bin/bash
# build Ipopt
set -euo pipefail
script_dir=$(dirname "$0")
script_dir=$(cd "$script_dir" && pwd )
cd "$script_dir"

make_opts=${MAKEOPTS-}

ipopt_url="https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.12.tgz"
ipopt_archive=Ipopt-3.12.12.tgz
ipopt_dir=$script_dir/Ipopt-3.12.12
ipopt_install_dir=$ipopt_dir/install

if [[ ! -d ${ipopt_dir} ]]; then
  if [[ ! -f ${ipopt_archive} ]]
    then wget -O $ipopt_archive ${ipopt_url}
  else
    echo "Skipping download step."
  fi

  tar -xzf ${ipopt_archive}

  # fetch dependencies
  cd "${ipopt_dir}"/ThirdParty/Blas
  ./get.Blas

  cd "${ipopt_dir}"/ThirdParty/Lapack
  ./get.Lapack

  cd "${ipopt_dir}"/ThirdParty/ASL
  ./get.ASL

  cd "${ipopt_dir}"/ThirdParty/HSL/

  if [[ ! -d coinhsl ]]; then
    # Need some coinhsl library. Check for common ones:
    if [[ -f "$script_dir"/coinhsl-2015.06.23.tar.gz ]]; then
      tar -xzf "$script_dir"/coinhsl-2015.06.23.tar.gz
      mv coinhsl-2015.06.23 coinhsl
    elif [[ -f "$script_dir"/coinhsl-2014.01.10.tar.gz ]]; then
                  tar -xzf "$script_dir"/coinhsl-2014.01.10.tar.gz
                  mv coinhsl-2014.01.10 coinhsl
    else
      echo "Did not find coinhsl/ or a known coinhsl archive."
      echo "Press any key to continue"
      read -n 1 -s -r
    fi
  fi

  cd "${ipopt_dir}"

  # With stupid case-insensitive macOS filesystem we cannot have INSTALL and install/
  if [[ -f install ]]; then
    mv INSTALL INSTALL.bak
  fi

  ./configure --prefix="$ipopt_install_dir" \
              --enable-static \
              --disable-shared \
              --with-pic
  make $make_opts
  make $make_opts install
else
    echo "Skipping building Ipopt. Directory already exists."
fi
