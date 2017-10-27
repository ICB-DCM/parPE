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
