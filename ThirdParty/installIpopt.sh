#!/bin/bash
# build Ipopt
set -e

tar -xzf Ipopt-3.12.9.tgz
cd Ipopt-3.12.9/

# fetch dependencies
cd ThirdParty/Blas
./get.Blas
cd ../Lapack
./get.Lapack
cd ../ASL
./get.ASL
cd ../HSL/
#TODO: need to get coinhsl-2015.06.23.tar.gz
tar -xzf ../../../coinhsl-2015.06.23.tar.gz
mv coinhsl-2015.06.23 coinhsl
cd ../..

./configure --prefix=`pwd`/install --enable-static --enable-shared
make -j 12
make install
