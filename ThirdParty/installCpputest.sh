#!/bin/bash
set -e 

tar -xzf cpputest-3.8.tar.gz
cd cpputest-3.8/
cd cpputest_build/
../configure --prefix=`pwd`/install
make -j12
make install
