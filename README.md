# Readme

## Building

Create a clean build folder

rm -rf build/; mkdir build ; cd build

To configure with MPI support rum cmake as

CC=mpicc CXX=mpic++ cmake ..

and build using 

make -j8

