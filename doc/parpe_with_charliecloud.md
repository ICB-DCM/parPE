# Using parPE with charliecloud

This document is an extension to [petab_model_import.md](petab_model_import.md),
describing how to perform parameter estimation using parPE inside a docker
container under charliecloud. 
[Charliecloud](https://hpc.github.io/charliecloud/) is a tool that "provides
user-defined software stacks (UDSS) for high-performance computing (HPC)
centers". We focus on charliecloud here, but it most content is easily
adaptable for use with plain docker.

## parPE base docker image

The next two sections describe how to build or download the default parPE
base docker image. The following section demonstrates how to perform parameter
estimation for an example model included in the the image.


### Generating parPE base docker image

This will create the parPE base image *from parPE from github*
(takes about 10'):

    cd container/charliecloud/parpe_base
    ch-build -t parpe

Export image to charliecloud archive in the current directory:

    ch-docker2tar parpe:latest .


### Using a provided docker image

Instead of building the docker image yourself, you can download a Ubuntu-based
parPE image from [dockerhub](https://hub.docker.com/r/dweindl/parpe) using:

    docker pull dweindl/parpe:latest


## Running the example model

This section show how to perform parameter estimation for an example model
included in the the image. The example here assume that you are using
interactive cluster access with SLURM. It should be easy to adapt the example
code to use a different submission system or to create a batch job.

### In an interactive session

Start an interactive session. With SLURM, e.g. 
`srun -p serial_fed28 --pty bash`

On the compute node:

    CHARLIE_DEST_DIR=/var/tmp
    CHARLIE_TAR=parpe\:latest.tar.gz
    OUTPUT_DIR=parpe_test # where results will be written to
    ch-tar2dir ${CHARLIE_TAR} ${CHARLIE_DEST_DIR}
    mkdir -p ${OUTPUT_DIR}
    ch-run -b ${OUTPUT_DIR}:/mnt/ ${CHARLIE_DEST_DIR}/parpe\:latest/ -- \
        mpirun /root/parPE/build/examples/parpeamici/steadystate/example_steadystate_multi \
        -o /mnt/pe-results/ /root/parPE/build/examples/parpeamici/steadystate/steadystate_scaled-prefix/src/steadystate_scaled/example_data.h5

The results will be written to `${OUTPUT_DIR}`.
