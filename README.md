# Readme

## Dependencies

* MPI
* PTHREADS
* IPOPT (coinhsl)
* CERES (Eigen)
* AMICI
* BLAS (cblas / intel mkl)

Scripts to fetch and build those third-party libraries are provided in `/ThirdParty/` 

## Building

Create a clean build folder

```
rm -rf build/; mkdir build ; cd build
```

To configure with MPI support rum cmake as
```
CC=mpicc CXX=mpiCC cmake ..
```
and build using 

```
make -j8
```

A sample build script is provided as `/install.sh`.


