#!/bin/bash
# build Ipopt
set -e
SCRIPT_DIR="`dirname \"$0\"`"
SCRIPT_DIR="`( cd \"$SCRIPT_DIR\" && pwd )`"
cd $SCRIPT_DIR

# unpack IpOpt source code
if [ ! -d Ipopt-3.12.9 ]; then
	tar -xzf Ipopt-3.12.9.tgz
fi

IPOPT_ROOT=$SCRIPT_DIR/Ipopt-3.12.9

# fetch dependencies
cd $IPOPT_ROOT/ThirdParty/Blas
./get.Blas

cd $IPOPT_ROOT/ThirdParty/Lapack
./get.Lapack

cd $IPOPT_ROOT/ThirdParty/ASL
./get.ASL

cd $IPOPT_ROOT/ThirdParty/HSL/

if [ ! -d coinhsl ]; then
	# Need some coinhsl library. Check for common ones:
	if [ -f ../../../coinhsl-2015.06.23.tar.gz ]; then
		tar -xzf ../../../coinhsl-2015.06.23.tar.gz
		mv coinhsl-2015.06.23 coinhsl
	elif [ -f ../../../coinhsl-2014.01.10.tar.gz ]; then
                tar -xzf ../../../coinhsl-2014.01.10.tar.gz
                mv coinhsl-2014.01.10 coinhsl
	else
		echo "Did not find coinhsl/ or a known coinhsl archive."
		echo "Press any key to continue"
		read -n 1 -s -r
	fi
fi

cd $IPOPT_ROOT

# With stupid case-insensitive macOS filesystem we cannot have INSTALL and install/
if [ -f install ]; then
	mv INSTALL INSTALL.bak
fi

./configure --prefix=`pwd`/install --enable-static --disable-shared --with-pic
make -j 12
make install
