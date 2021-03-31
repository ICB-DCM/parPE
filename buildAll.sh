#!/bin/bash

set -euo pipefail

parpe_root=$(dirname "$0")
parpe_root=$(cd "${parpe_root}" && pwd)

make_opts=${MAKEOPTS--j12}

# Build amici
amici_dir="${parpe_root}/deps/AMICI/"
cd "${amici_dir}"
#scripts/buildAll.sh

# build ceres
"${parpe_root}/ThirdParty/installCeres.sh"
# build ipopt
"${parpe_root}/ThirdParty/installIpopt.sh"

# build parpe
export PKG_CONFIG_PATH=${parpe_root}/ThirdParty/Ipopt-releases-3.13.3/:${PKG_CONFIG_PATH:-}
parpe_build_dir="${parpe_root}/build"
mkdir -p "${parpe_build_dir}"
cd "${parpe_build_dir}"
cmake "${parpe_root}"
make ${make_opts}

# run tests
make test
