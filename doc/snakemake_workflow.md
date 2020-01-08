# Using parPE with Snakemake

*parPE* comes with some rudimentary 
[Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows.


## Parameter estimation

A very basic Snakefile for parameter estimation is located in `snakemake/`.
This workflow is adapted to a specific model via a YAML configuration file.
An example is shown in 
[`snakemake/parpe_optimize_petab_steadystate.yaml`](../snakemake/parpe_optimize_petab_steadystate.yaml).


### An example

If parPE was build with examples, this workflow can be run after adapting
some paths in the example config file. `petab.root` needs to be adapted to
reflect your parPE build directory, and probably you want to change the output
directories as well. The documented schema for this configuration file is
provided in `snakemake/config.schema.yaml`.

After that you can run the full pipeline with:

    cd snakemake
    snakemake --configfile parpe_optimize_petab_steadystate.yaml -- postprocess

This generate C++ code of the model, build model specific binaries for
parameter estimation, run parameters, and process the results.
After successful completion, you should see a file `test.png` in your working
directory showing the optimization trajectory (as well as other files, once
this gets extended).

The different snakemake targets are explained in the Snakefile and can be
accessed using `snakemake --list` from the directory in which the Snakefile
is located.

To use this workflow with a different PEtab problem, you should only need to
adapt the configuration file. 
