#!/usr/bin/env bash
# build Ipopt

set -euo pipefail
script_dir=$(dirname "$0")
script_dir=$(cd "${script_dir}" && pwd )
cd "${script_dir}"

make_opts=${MAKEOPTS-}

git clone https://github.com/coin-or-tools/ThirdParty-HSL.git
cd ThirdParty-HSL
# Now unpack the HSL sources archive, move and rename the resulting directory so that it becomes ThirdParty-HSL/coinhsl. Then, in ThirdParty-HSL, configure, build, and install the HSL sources:
./configure
make
sudo make install


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
