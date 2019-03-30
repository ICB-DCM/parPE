#!/bin/bash
# Download and build googletest

set -e

SCRIPT_PATH="`dirname \"$0\"`"
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"

cd ${SCRIPT_PATH}

if [[ ! -d googletest ]]; then
  echo "googletest/ does not exist"
  if [[ ! -f googletest.zip ]]; then
    echo "Downloading googletest..."
    wget https://github.com/google/googletest/archive/master.zip -O googletest.zip
  fi

  echo "Unpacking and building googletest..."
  unzip googletest
  mv googletest-master googletest
  cd googletest
  mkdir -p build
  cd build
  cmake ..
  make -j 4
else
  echo "googletest/ exists. nothing to do."
fi
