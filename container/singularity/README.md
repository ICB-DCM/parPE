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


TO BE EXTENDED
