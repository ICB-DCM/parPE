# Running parameter optimization for an AMICI model using `parpe::OptimizationApplication`

This document describes the output and options of parPE executables based on
`parpe::OptimizationApplication`.

**TODO** 

For now see example notebook in `../examples/parpeamici/steadystate/`.

## Command line options

Run the created executable with the `-h`/`--help` argument.

## Output explained

### Optimizer iterations 

```
[2018-08-26 12:23:02] [INF] [   0/i21r01c01s01] [o4i2] iter: 2 cost: 50339.6 time_iter: wall: 278.65s cpu: 59638.5s time_optim: wall: 1312.94s cpu: 230120s
 ^timestamp            ^loglvl ^rank/hostname ^optimizerProgress
```

- Log level: one of debug ('DBG'), information ('INF'), warning ('WRN'), error ('ERR'), critical ('CRI')
- Rank: MPI rank. Master is 0. For non-MPI runs this is -1. 
- Optimizer progress: Number after 'o' a running index of optimizer runs, number after 'i' is the current iteration
- 'iter': current iteration
- 'cost': current cost function value
- 'time_iter: wall:' wall-time for this iteration, 'cpu' approximate CPU-time for this iteration
- 'time_optim': cumulative 'time_iter' for the given optimizer run

### Simulations during cost function evaluation  

```
[2018-09-05 12:14:59] [DBG] [18/i20r01c01s11] [o0i0c1717] Result for 2379: -0,792136 (0) (21,6578s+)
 ^timestamp           ^loglvl ^rank/hostname ^optimizerProgress   ^messageID ^llh ^status ^time ^gradientFlag
```
- *see above*
- Optimizer progress: *see above* + Number after 'c' is the current condition index
- messageID: internal use only
- llh: Log-likelihood as provided by AMICI
- status: AMICI status
- time: simulation wall-time
- gradientFlag: '-' indicates simulation without gradient computation;
  'A' indicates gradient computation using adjoint sensitivity analysis;
  'F' indicates gradient computation using forward sensitivity analysis;
  this is followed by a number indicating w.r.t. how many parameters the
  gradient was computed.

## Environment variables

- **PARPE_NUM_PARALLEL_STARTS**

  Setting `PARPE_NUM_PARALLEL_STARTS=n` will create a maximum of `n` threads
  for concurrent local optimizations during multi-start optimization.
  If unset, this defaults to the number of concurrent threads supported by
  hardware.

- **PARPE_LOG_SIMULATIONS**

  Setting `PARPE_LOG_SIMULATIONS=1` will cause every single AMICI simulation to be saved in the result files.
  
  Note: In the case of large models and many simulation conditions this can result in huge files. 

- **PARPE_NO_DEBUG**

  With `PARPE_NO_DEBUG=1` no `LOGLVL_DEBUG` messages (i.e. those prefixed with `[DBG]`) will be printed.  

- **PARPE_MAX_SIMULATIONS_PER_PACKAGE=n** and **PARPE_MAX_GRADIENT_SIMULATIONS_PER_PACKAGE=n**
  
  These environment variables determine how many simulation requests will be sent to workers within one message.
  For simulations without sensitivities (e.g. during IpOpt line-search) `PARPE_MAX_SIMULATIONS_PER_PACKAGE` is used,
  for simulations with sensitivities `PARPE_MAX_GRADIENT_SIMULATIONS_PER_PACKAGE` is used.
  Changing these values may be useful to reduce communication overhead when working with a high number of workers.
  It may also be necessary to reduce this value in the case only very few simulations can be performed in parallel.
  
  Note: These variables have no effect in case of shared-memory (non-MPI) execution
   
- **PARPE_NUM_SIMULATION_TRIALS** (integer) and
  **PARPE_INTEGRATION_TOLERANCE_RELAXATION_FACTOR** (float)

  In case of simulation failure, parPE try rerunning an AMICI simulation with a
  `PARPE_INTEGRATION_TOLERANCE_RELAXATION_FACTOR`-fold higher error tolerance
  for a total of `PARPE_NUM_SIMULATION_TRIALS` times (including the initial
  attempt).
