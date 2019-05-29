[![Run Status](https://api.shippable.com/projects/59463d3e8993d7070010407b/badge?branch=master)](https://app.shippable.com/github/dweindl/parPE)
[![Coverage Badge](https://api.shippable.com/projects/59463d3e8993d7070010407b/coverageBadge?branch=master)](https://app.shippable.com/github/dweindl/parPE)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/1f1ee5a0d90d431499f200a148fb7fdc)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ICB-DCM/parPE&amp;utm_campaign=Badge_Grade)
# parPE

The *parPE* library provides functionality for solving large-scale parameter
optimization problems requiring up to thousands of simulations per objective
function evaluation on HPC systems.

Currently, parPE provides interfaces to the
[Ipopt](http://www.coin-or.org/Ipopt/),
[Ceres](http://ceres-solver.org/),
[FFSQP](https://www.isr.umd.edu/news/news_story.php?id=4088) and
[SUMSL (CALGO/TOMS 611)](http://www.netlib.org/toms/index.html)
optimizers. parPE offers easy integration with
[AMICI](https://github.com/ICB-DCM/AMICI)-generated ordinary differential
equation (ODE) models.

## Features

parPE offers the following features:

* MPI-based load-balancing of individual simulations
* improved load balancing by intermingling multiple optimization runs
  (multi-start local optimization)
* simple integration with [SBML](http://sbml.org/) models via
  [AMICI](https://github.com/ICB-DCM/AMICI) and
  [PEtab](https://github.com/ICB-DCM/PEtab)
* interfaces to Ipopt and Ceres optimizers
* HDF5 I/O compatible with a wide variety of programming languages
* Good parallel scaling to up to several thousand cores
  (highly problem dependent)

## Dependencies

For full functionality, parPE requires the following libraries:

* CMAKE (>=3.6)
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
* C++14 compiler
* Python >= 3.6, including header files

On Debian-based systems, dependencies can be installed via:
```
sudo apt-get install build-essential gfortran libmpich-dev libblas-dev \
    libhdf5-dev cmake libceres-dev coinor-libipopt-dev libcpputest-dev \
    libboost-serialization-dev libpython-dev
```

Scripts to fetch and build the remaining dependencies are provided in
`/ThirdParty/`:

```
cd ThirdParty
./downloadPackages.sh
./installDeps.sh
```

NOTE: When using `ThirdParty/installIpopt.sh` to build Ipopt, follow the
instructions in `ThirdParty/Ipopt-3.12.12/ThirdParty/HSL/INSTALL.HSL` for
obtaining the hsl library before continuing, otherwise IpOpt will not be
usable. Afterwards, (re)run `ThirdParty/installIpopt.sh`.


## Building

After having taken care of the dependencies listed above, parPE can be built:

```
./build.sh
```

Other sample build scripts are provided as `/build*.sh`.

## Recently tested compilers

* GCC 8.3.0
* Intel icpc (ICC) 17.0.6


## Documentation & further information

No extensive full-text documentation is available yet. See `*/examples` and
`*/tests` for usage examples. Some additional documentation is available
in `doc/` and among [Github issues](https://github.com/ICB-DCM/parPE/issues).

## FAQ

Q: The program is killed due to memory exhaustion, what should I do?

A: When running with MPI, the master process (rank 0) is consuming more memory
than the others. Consider reserving more memory for this one. For LoadLeveler,
this can be done conveniently via `#@ first_node_tasks`.

## References

parPE has been used in the following projects:

- Leonard Schmiester, Yannik Schälte, Fabian Fröhlich, Jan Hasenauer, Daniel Weindl.
  *Efficient parameterization of large-scale dynamic models based on relative measurements*.
  bioRxiv 2018.
  [doi:10.1101/579045](https://www.biorxiv.org/content/10.1101/579045v1)

- [CanPathPro](http://canpathpro.eu/)
