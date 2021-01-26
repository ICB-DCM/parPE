#!/usr/bin/env bash
# build Ipopt

set -euo pipefail
script_dir=$(dirname "$0")
script_dir=$(cd "${script_dir}" && pwd )
cd "${script_dir}"

make_opts=${MAKEOPTS-}

ipopt_url="https://www.coin-or.org/download/source/Ipopt/Ipopt-3.13.3.tgz"
ipopt_archive="Ipopt-3.13.3.tgz"
ipopt_dir="${script_dir}/Ipopt-releases-3.13.3"
ipopt_install_dir="${ipopt_dir}/install"

if [[ ! -d "${ipopt_dir}" ]]; then
  if [[ ! -f "${ipopt_archive}" ]]
    then wget -O "${ipopt_archive}" "${ipopt_url}"
  else
    echo "Skipping download step."
  fi

  tar -xzf "${ipopt_archive}"

  cd "${ipopt_dir}"

  ./configure --prefix="${ipopt_install_dir}" \
              --enable-static \
              --disable-shared \
              --with-pic
  make $make_opts
  make $make_opts install
else
    echo "Skipping building Ipopt. Directory already exists."
fi

# Handle HSL solvers
cd "${ipopt_dir}"
if [[ ! -d "ThirdParty-HSL" ]]; then
  git clone https://github.com/coin-or-tools/ThirdParty-HSL.git
fi

cd ThirdParty-HSL

if [[ ! -d "coinhsl" ]]; then
  # Need some coinhsl library. Check for common ones:
  if [[ -f "${script_dir}/coinhsl-2015.06.23.tar.gz" ]]; then
    tar -xzf "${script_dir}/coinhsl-2015.06.23.tar.gz"
    mv coinhsl-2015.06.23 coinhsl
  elif [[ -f "${script_dir}/coinhsl-2014.01.10.tar.gz" ]]; then
                tar -xzf "${script_dir}/coinhsl-2014.01.10.tar.gz"
                mv coinhsl-2014.01.10 coinhsl
  else
    echo "Did not find coinhsl/ or a known coinhsl archive."
    echo "Press any key to continue"
    read -n 1 -s -r
  fi
fi

./configure --prefix="$(pwd)/install"
make
make install
(cd install/lib && test ! -e libhsl.so && ln -s libcoinhsl.so libhsl.so || true)

