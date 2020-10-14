#!/usr/bin/env bash
# Install CBLAS

set -euo pipefail
script_dir=$(dirname "$0")
script_dir=$(cd "$script_dir" && pwd )

make_opts=${MAKEOPTS-}
cd "$script_dir"

cblas_url="http://www.netlib.org/blas/blast-forum/cblas.tgz"
cblas_archive="cblas.tgz"
cblas_dir="${script_dir}/CBLAS"

if [[ ! -d "${cblas_dir}" ]]; then
    if [[ ! -f "${cblas_archive}" ]]
        then wget -O "${cblas_archive}" "${cblas_url}"
    fi

    tar -xzf "${cblas_archive}"
fi

cd "$cblas_dir"
ln -sf Makefile.LINUX Makefile.in
make $make_opts alllib

