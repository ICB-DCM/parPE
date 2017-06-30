 [![Run Status](https://api.shippable.com/projects/59463d3e8993d7070010407b/badge?branch=master)](https://app.shippable.com/github/dweindl/parPE)
 [![Coverage Badge](https://api.shippable.com/projects/59463d3e8993d7070010407b/coverageBadge?branch=master)](https://app.shippable.com/github/dweindl/parPE) 

# parPE

The *parPE* library provides functionality for solving large-scale parameter optimization
problems requiring thousands of simulations per objective function evaluation on HPC systems.

Currently, parPE provides interfaces to IpOpt and Ceres optimizers. parPE offers easy integration with
AMICI-generated ordinary differential equation (ODE) models.

## Features

parPE offers the following features:

* MPI-based load-balancing of individual simulations
* improved load balancing by intermingeling multiple optimization runs (multi-start local optimization)
* simple integration with SBML models via AMICI
* interfaces to Ipopt and Ceres optimizers
* HDF5 I/O compatible with a wide variety of programming languages
* Good parallel scaling to up to several thousand cores (highly problem dependent)

## Dependencies

For full functionality, parPE requires the following libraries:

* MPI
* PTHREADS
* IPOPT (requires coinhsl)
* CERES (requires Eigen)
* HDF5
* BLAS (CBLAS / Intel MKL)
* AMICI [https://github.com/ICB-DCM/AMICI](AMICI) (requires SuiteSparse, Sundials)

Scripts to fetch and build those third-party libraries are provided in `/ThirdParty/` 

On Debian-based systems, dependencies can be installed via:
```
sudo apt-get install build-essential gfortran libmpich-dev libblas-dev libhdf5-dev cmake \
    libceres-dev coinor-libipopt-dev libcpputest-dev
```

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

Sample build scripts are provided as `/build*.sh`.

## Tested compilers

* GCC 4.8.3
* GCC 5.4.0
* Intel ...


## Documentation & further information

No extensive full-text documentation is available yet. See `*/examples` and `*/tests` for usage examples. 
Little documentation is available in `doc` and among github issues. 
