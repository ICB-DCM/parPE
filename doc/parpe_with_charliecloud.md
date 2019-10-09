# Using parPE with charliecloud

This document is an extension to [petab_model_import.md](petab_model_import.md),
describing how to perform parameter estimation using parPE inside a docker
container under charliecloud. 
[Charliecloud](https://hpc.github.io/charliecloud/) is a tool that "provides
user-defined software stacks (UDSS) for high-performance computing (HPC)
centers". We focus on charliecloud here, but most of the steps are easily
adapted for use with plain docker.

## parPE base docker image

The next two sections describe how to build or download the default parPE
base docker image. The following section demonstrates how to perform parameter
estimation for an example model included in the the image.


### Generating parPE base docker image

This will create the parPE base image *from parPE from github*
(takes about 10'):

    cd container/charliecloud/parpe_base
    ch-build -t parpe .

Export image to charliecloud archive in the current directory:

    ch-docker2tar parpe:latest .


### Using a provided docker image

Instead of building the docker image yourself, you can download a Ubuntu-based
parPE image from [dockerhub](https://hub.docker.com/r/dweindl/parpe) using:

    docker pull dweindl/parpe:latest

**Developer note:**
To update the image on dockerhub, run:

    sudo docker images # check for IMAGE ID
    sudo docker tag $IMAGE_ID dweindl/parpe:latest
    sudo docker push dweindl/parpe:latest 

(TODO: auto-deploy after CI)

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
    ch-tar2dir "${CHARLIE_TAR}" "${CHARLIE_DEST_DIR}"
    mkdir -p "${OUTPUT_DIR}"
    ch-run -b "${OUTPUT_DIR}":/mnt/ "${CHARLIE_DEST_DIR}"/parpe\:latest/ -- \
        mpirun /root/parPE/build/examples/parpeamici/steadystate/example_steadystate_multi \
        -o /mnt/pe-results/ \
        /root/parPE/build/examples/parpeamici/steadystate/steadystate_scaled-prefix/src/steadystate_scaled/example_data.h5

The results will be written to `${OUTPUT_DIR}`.


## Optimizing a PEtab model

1. Generate / fetch a parPE docker image as described above

1. Generate charliecloud image (potentially requires `sudo`)
 
    `ch-docker2tar parpe:latest .`

1. Extract charliecloud image (adapt paths as needed)

    ```
    CHARLIE_DEST_DIR=/var/tmp
    CHARLIE_TAR=parpe\:latest.tar.gz
    ch-tar2dir "${CHARLIE_TAR}" "${CHARLIE_DEST_DIR}"
    ```

1. Place your PEtab files into `$PETAB_DIR`

   For testing you can download an example model via:

    ```    
    wget "https://raw.githubusercontent.com/LeonardSchmiester/Benchmark-Models/hackathon/hackathon_contributions_new_data_format/Zheng_PNAS2012/model_Zheng_PNAS2012.xml"
    wget "https://raw.githubusercontent.com/LeonardSchmiester/Benchmark-Models/hackathon/hackathon_contributions_new_data_format/Zheng_PNAS2012/measurementData_Zheng_PNAS2012.tsv"
    wget "https://raw.githubusercontent.com/LeonardSchmiester/Benchmark-Models/hackathon/hackathon_contributions_new_data_format/Zheng_PNAS2012/experimentalCondition_Zheng_PNAS2012.tsv"
    wget "https://raw.githubusercontent.com/LeonardSchmiester/Benchmark-Models/hackathon/hackathon_contributions_new_data_format/Zheng_PNAS2012/parameters_Zheng_PNAS2012.tsv"
    ```

1. Create Snakemake YAML config for those files

    For the example files above, run:
    
    ```
    cat > parpe_optimize_petab.yaml << EOF
   petab:
        sbml_file: 'model_Zheng_PNAS2012.xml'
        measurement_file: 'measurementData_Zheng_PNAS2012.tsv'
        condition_file: 'experimentalCondition_Zheng_PNAS2012.tsv'
        parameter_file: 'parameters_Zheng_PNAS2012.tsv'
    
    model_name: 'Zheng_PNAS2012'
    amici_build_dir: '/root/parPE/deps/AMICI/build'
    amici_src_dir: '/root/parPE/deps/AMICI/'
    parpe_build_dir: '/root/parPE/build/'
    EOF
    ```

1. Run the Snakemake workflow, e.g. as

    ```    
    PETAB_DIR=Zheng_PNAS2012 # where results will be written to
    SNAKEMAKE_CONFIG=/mnt/parpe_optimize_petab.yaml
    mkdir -p "${PETAB_DIR}"
    ch-run -b "${PETAB_DIR}":/mnt/ --no-home -c /mnt/ \
        --unset-env=AMICI_ROOT \
        "${CHARLIE_DEST_DIR}"/parpe\:latest/ -- \
        /root/parPE/misc/run_in_venv.sh /root/parPE/build/venv snakemake \
        -s /root/parPE/snakemake/Snakefile \
        --configfile "${SNAKEMAKE_CONFIG}" -- postprocess
    ```

    For more details on this workflow, see
    [../snakemake/Snakefile](../snakemake/Snakefile).
