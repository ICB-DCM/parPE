# build EIGEN
tar -xzf eigen-3.3.3.tar.gz
cd eigen-eigen-67e894c6cd8f
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=`pwd`/install ..
make -j12
make install
cd ../..

# build CERES
tar -xzf ceres-solver-1.12.0.tar.gz
cd ceres-solver-1.12.0/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=`pwd`/install \
      -DBUILD_SHARED_LIBS=ON \
      -DBUILD_TESTING=OFF \
      -DGFLAGS=OFF \
      -DEIGEN_INCLUDE_DIR=`pwd`/../../eigen-eigen-67e894c6cd8f/ \
      -DMINIGLOG=ON \
      ..
      
make -j12
make install

      
