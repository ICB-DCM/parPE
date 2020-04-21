#!/usr/bin/env bash
# build f2c

set -euo pipefail

f2c_archive=libf2c.zip
if [ ! -d libf2c ]; then
  if [ ! -f "${f2c_archive}" ]; then
    wget "http://www.netlib.org/f2c/libf2c.zip" -o "${f2c_archive}"
  fi
    unzip "${f2c_archive}" -d libf2c
fi

cd libf2c
cp makefile.u Makefile

make hadd
make -j12

cd ../toms611/
ln -sf ../libf2c/f2c.h .
