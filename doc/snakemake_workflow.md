# Using parPE with Snakemake

*parPE* comes with some rudimentary 
[Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow.


## Parameter estimation

A very basic Snakefile for parameter estimation is located in `snakemake/`.
This workflow is adapted to a specific model via a YAML configuration file.
An example is shown in 
[`snakemake/parpe_optimize_petab_steadystate.yaml`](../snakemake/parpe_optimize_petab_steadystate.yaml).


### An example

If parPE was built with examples, this workflow can be run after adapting
some paths in the example config file. `petab.root` needs to be adapted to
reflect your parPE build directory, and probably you want to change the output
directories as well. The documented schema for this configuration file is
provided in `snakemake/config.schema.yaml`.

After that you can run the full workflow with:

```shell
cd snakemake
snakemake --configfile parpe_optimize_petab_steadystate.yaml -- postprocess
```

This generates C++ code of the model, builds model specific binaries for
parameter estimation, runs parameter estimation, and processes the results.
After successful completion, you should see a file
`results/figures/cost_trajectory.png` showing the optimizer trajectory.

The different snakemake targets are explained in the Snakefile and can be
accessed using `snakemake --list` from the directory in which the `Snakefile`
is located.

To use this workflow with a different PEtab problem, you should only need to
adapt the configuration file. 
