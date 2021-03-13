<a href="https://github.com/ICB-DCM/parPE/actions?query=workflow%3A%22parPE+tests%22">
<img src="https://github.com/ICB-DCM/parPE/workflows/parPE%20tests/badge.svg?branch=master" alt="parPE tests"></a>
<a href="https://sonarcloud.io/dashboard?id=ICB-DCM_parPE">
<img src="https://sonarcloud.io/api/project_badges/measure?project=ICB-DCM_parPE&metric=coverage" alt="Coverage"></a>
<a href="https://github.com/ICB-DCM/parPE/actions?query=workflow%3A%22PEtab+test+suite%22">
<img src="https://github.com/ICB-DCM/parPE/workflows/PEtab%20test%20suite/badge.svg" alt="PEtab test suite"></a>
<a href="https://github.com/ICB-DCM/parPE/actions?query=workflow%3A%22Deploy+to+dockerhub%22">
<img src="https://github.com/ICB-DCM/parPE/workflows/Deploy%20to%20dockerhub/badge.svg" alt="Deploy to dockerhub"></a>
<a href="https://doi.org/10.5281/zenodo.3478612">
<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3478612.svg" alt="DOI"></a>

# parPE

The *parPE* library provides functionality for solving large-scale parameter
optimization problems requiring up to thousands of simulations per objective
function evaluation on high performance computing (HPC) systems.

parPE offers easy integration with
[AMICI](https://github.com/ICB-DCM/AMICI)-generated ordinary differential
equation (ODE) models.

## Features

parPE offers the following features:

* MPI-based load-balancing of individual simulations
* improved load balancing by intermingling multiple optimization runs
  (multi-start local optimization)
* simple integration with [SBML](http://sbml.org/) models via
  [AMICI](https://github.com/ICB-DCM/AMICI) and
  [PEtab](https://github.com/PEtab-dev/PEtab)
* interfaces to [Ipopt](http://www.coin-or.org/Ipopt/),
  [Ceres](http://ceres-solver.org/),
  [FFSQP](https://www.isr.umd.edu/news/news_story.php?id=4088) and
  [SUMSL (CALGO/TOMS 611)](http://www.netlib.org/toms/index.html) optimizers
* HDF5 I/O compatible with a wide variety of programming languages
* Good parallel scaling to up to several thousand cores
  (highly problem dependent)

## Getting started

Although various modules of parPE can be used independently, the most
meaningful and convenient use case is parameter optimization for an SBML model
specified in the [PEtab](https://github.com/ICB-DCM/PEtab) format. This is
described in [doc/petab_model_import.md](doc/petab_model_import.md).

## Dependencies

For full functionality, parPE requires the following libraries:

* CMAKE (>=3.10)
* MPI ([OpenMPI](https://www.open-mpi.org/),
  [MPICH](https://www.mpich.org/), ...)
* PTHREADS
* IPOPT (>= 1.2.7) (requires coinhsl)
* CERES (>=1.13)
  ([requires Eigen](http://ceres-solver.org/installation.html#dependencies))
* HDF5 (>= 1.10)
* CBLAS compatible BLAS (libcblas, Intel MKL, ...)
* [AMICI](https://github.com/ICB-DCM/AMICI) (included in this repository)
  (uses SuiteSparse, Sundials)
* C++17 compiler
* Python >= 3.7, including header files

On Debian-based systems, dependencies can be installed via:
```shell
sudo apt-get install build-essential cmake cmake-curses-gui \
    coinor-libipopt-dev curl gfortran \
    libblas-dev libboost-serialization-dev libceres-dev libcpputest-dev \
    libmpich-dev libhdf5-dev libpython3-dev python3-pip
```

Scripts to fetch and build the remaining dependencies are provided in
`/ThirdParty/`:

```shell
ThirdParty/installDeps.sh
```

NOTE: When using `ThirdParty/installIpopt.sh` to build Ipopt, you may have to
download the HSL library separately as described at
https://coin-or.github.io/Ipopt/INSTALL.html#DOWNLOAD_HSL. Place the HSL
archive into `ThirdParty` before running `ThirdParty/installIpopt.sh`. If asked
type in your coinhsl version (e.g. `2019.05.21` if you
have `coinhsl-2019.05.21.tar.gz`).


## Building

After having taken care of the dependencies listed above, parPE can be built:

```shell
./buildAll.sh
```

Other sample build scripts are provided as `/build*.sh`.

## Recently tested compilers

* GCC 10.2.0
* Intel icpc (ICC) 17.0.6

## Docker

There is a Dockerfile available in `container/charliecloud/` and images
can be found on [dockerhub](https://hub.docker.com/r/dweindl/parpe/).

## Documentation & further information

Some high-level documentation is available at
[https://parpe.readthedocs.io/en/latest/](https://parpe.readthedocs.io/en/latest/)
and among [Github issues](https://github.com/ICB-DCM/parPE/issues). No extensive
full-text documentation is available for the C++ interface yet. For usage of
the C++ interface see [`examples/`](examples/) and `*/tests`.


## References

parPE is being used or has been used in the following projects:

- Leonard Schmiester, Yannik Schälte, Fabian Fröhlich, Jan Hasenauer,
  Daniel Weindl.
  *Efficient parameterization of large-scale dynamic models based on relative measurements*.
  Bioinformatics, btz581, [doi:10.1093/bioinformatics/btz581](https://doi.org/10.1093/bioinformatics/btz581)
  (preprint: [doi:10.1101/579045](https://www.biorxiv.org/content/10.1101/579045v1)).

- Paul Stapor, Leonard Schmiester, Christoph Wierling, Bodo Lange,
  Daniel Weindl, and Jan Hasenauer. 2019.
  *Mini-Batch Optimization Enables Training of Ode Models on Large-Scale Datasets.*
  bioRxiv. Cold Spring Harbor Laboratory.
  preprint: [doi:10.1101/859884](https://doi.org/10.1101/859884).

- [CanPathPro](http://canpathpro.eu/)


## Funding

parPE has been developed within research projects receiving external funding:

* Through  the  European  Union's  Horizon  2020  research  and innovation
  programme under grant agreement no. 686282
  ([CanPathPro](http://canpathpro.eu/)).

* Computer resources for testing parPE have been provided among others by the 
  Gauss Centre for Supercomputing / Leibniz Supercomputing Centre under grant
  pr62li.
