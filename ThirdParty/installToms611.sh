#!/bin/bash

set -e 

# build f2c

F2C_ARCHIVE=libf2c.zip
if [ ! -d libf2c ]; then
  if [ ! -f $F2C_ARCHIVE ]; then
    wget "http://www.netlib.org/f2c/libf2c.zip" -o $DLIB_ARCHIVE
  fi
    unzip $F2C_ARCHIVE -d libf2c
fi 

cd libf2c
cp makefile.u Makefile

make hadd
make -j12

cd ../toms611/
ln -sf ../libf2c/f2c.h
