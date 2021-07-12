#!/usr/bin/env bash
# Download and build DLIB

set -euo pipefail
script_dir=$(dirname "$0")
script_dir=$(cd "${script_dir}" && pwd)
cd "${script_dir}"

make_opts=${MAKEOPTS-}

dlib_archive="dlib-19.7.tar.bz2"
dlib_dir=${script_dir}/dlib-19.7

if [ ! -d "${dlib_dir}" ]; then
  if [ ! -f "${dlib_archive}" ]; then
    wget "http://dlib.net/files/${dlib_archive}" -O "${dlib_archive}"
  fi
  tar -xjf "${dlib_archive}"
fi

cd "${dlib_dir}"
mkdir -p build && cd build/
cmake -S .. \
  -DCMAKE_INSTALL_PREFIX="$(pwd)/install" \
  -DDLIB_NO_GUI_SUPPORT=ON \
  -DDLIB_GIF_SUPPORT=OFF \
  -DDLIB_JPEG_SUPPORT=OFF \
  -DDLIB_PNG_SUPPORT=OFF \
  -DDLIB_LINK_WITH_SQLITE3=OFF \
  ..

make ${make_opts} && make install
