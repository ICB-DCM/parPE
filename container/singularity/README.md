# Using parPE with singularity

## Singularity

> Singularity is an open source container platform designed to be simple, fast,
> and secure. Singularity is optimized for compute focused enterprise and HPC
>workloads, allowing untrusted users to run untrusted containers in a trusted
>way.

--- [https://github.com/hpcng/singularity](https://github.com/hpcng/singularity)

Documentation: [https://sylabs.io/guides/3.0/user-guide/index.html](https://sylabs.io/guides/3.0/user-guide/index.html)


## Using parPE with singularity

Singularity images can be created from available docker containers using:

```
singularity pull docker://dweindl/parpe:develop
```

To create a custom docker containers, see 
https://parpe.readthedocs.io/en/latest/parpe_with_charliecloud.html#generating-parpe-base-docker-image


An example for parameter estimation for the model at https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab/tree/master/Benchmark-Models/Zheng_PNAS2012 :

```shell
AMICI_MODEL_DIR=Zheng_PNAS2012
PARPE_MODEL_DIR=Zheng_PNAS2012_parpe
PETAB_YAML_FILE=Benchmark-Models-PEtab/Benchmark-Models/Zheng_PNAS2012/Zheng_PNAS2012.yaml
MODEL_NAME=Zheng_PNAS2012
H5_PE_INPUT=$PARPE_MODEL_DIR/data.h5

singularity exec parpe_develop.sif /parPE/misc/run_in_venv.sh /parPE/build/venv amici_import_petab -v -y $PETAB_YAML_FILE &
singularity exec parpe_develop.sif /parPE/misc/setup_amici_model.sh ${AMICI_MODEL_DIR} ${PARPE_MODEL_DIR}
singularity exec parpe_develop.sif /parPE/misc/run_in_venv.sh /parPE/build/venv parpe_petab_to_hdf5 \
    -n ${MODEL_NAME} \
    -y ${PETAB_YAML_FILE} \
    -d ${AMICI_MODEL_DIR} \
    -o ${H5_PE_INPUT}
 singularity exec parpe_develop.sif ${PARPE_MODEL_DIR}/build/estimate_$MODEL_NAME -o ${MODEL_NAME}_results/ $H5_PE_INPUT
```
