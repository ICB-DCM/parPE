 [![Run Status](https://api.shippable.com/projects/59463d3e8993d7070010407b/badge?branch=master)](https://app.shippable.com/github/dweindl/parPE)
 [![Coverage Badge](https://api.shippable.com/projects/59463d3e8993d7070010407b/coverageBadge?branch=master)](https://app.shippable.com/github/dweindl/parPE) 

# parPE

The *parPE* library provides functionality for solving large-scale parameter optimization
problems requiring thousands of simulations per objective function evaluation on HPC systems.

Currently, parPE provides interfaces to [IpOpt](http://www.coin-or.org/Ipopt/) and [Ceres](http://ceres-solver.org/) optimizers. parPE offers easy integration with
[AMICI](https://github.com/ICB-DCM/AMICI)-generated ordinary differential equation (ODE) models.

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

* CMAKE
* MPI
* PTHREADS
* IPOPT (requires coinhsl)
* CERES (requires Eigen)
* HDF5 (>= 1.10)
* BLAS (CBLAS / Intel MKL)
* [AMICI](https://github.com/ICB-DCM/AMICI) (requires SuiteSparse, Sundials)

On Debian-based systems, dependencies can be installed via:
```
sudo apt-get install build-essential gfortran libmpich-dev libblas-dev libhdf5-dev cmake \
    libceres-dev coinor-libipopt-dev libcpputest-dev libboost-serialization-dev
```

Scripts to fetch and build the remaining dependencies are provided in `/ThirdParty/` :

```
cd ThirdParty
./downloadPackages.sh
./installDeps.sh
```

## Building

After having taken care of the dependencies listed above, parPE can be built: 

```
./build.sh
```

Other sample build scripts are provided as `/build*.sh`.

## Tested compilers

* GCC 4.8.3
* GCC 5.4.0
* Intel 16.0.4


## Documentation & further information

No extensive full-text documentation is available yet. See `*/examples` and `*/tests` for usage examples. 
Little documentation is available in `doc` and among github issues. 
