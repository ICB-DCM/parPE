#!/bin/bash

set -e 

# build DLIB
DLIB_ARCHIVE="dlib-19.7.tar.bz2"
if [ ! -d dlib-19.7 ]; then
  if [ ! -f $DLIB_ARCHIVE ]; then
    wget "http://dlib.net/files/$DLIB_ARCHIVE" -o $DLIB_ARCHIVE
  fi
  tar -xjf $DLIB_ARCHIVE
fi 

cd dlib-19.7/
mkdir -p build && cd build/
 cmake -DCMAKE_INSTALL_PREFIX="`pwd`/install" \
    -DDLIB_NO_GUI_SUPPORT=ON \
    -DDLIB_GIF_SUPPORT=OFF \
    -DDLIB_JPEG_SUPPORT=OFF \
    -DDLIB_PNG_SUPPORT=OFF \
    -DDLIB_LINK_WITH_SQLITE3=OFF \
    ..
 
make -j12 && make install
