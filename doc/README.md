# parPE Documentation

This directory contains the parPE documentation.

Contents in order of assumed relevance:

- petab_model_import.md: Guide for using parPE for parameter estimation of a
  PEtab problem

- hdf5.md: Description of HDF5 file formats

- optimizationApplication.md: Documentation stub for using parameter estimation
  executables based on class `parpe::OptimizationApplication`, describing
  invocation, environment variables, and output

- standaloneSimulator.md: Documentation for simulator executable

- FAQ.md: Frequently asked questions

- Doxyfile.in: Configuration file for doxygen C++ documentation generation

- gfx/: Graphics for the documentation


# Building documentation

To build doxygen documentation run `doxygen` inside the `doc` directory.

To build the sphinx documentation after having built the doxygen documentation,
run `make html` from inside the `build/venv/` python environment. 

View the result with `firefox build/html/index.html`.
