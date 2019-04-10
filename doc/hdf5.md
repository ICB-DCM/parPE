# HDF5 data structure

```
+ /amiciOptions/
+ /optimizationOptions/
+ /fixedParameters/
+ /measurements/
+ /parameters/
```

## /amiciOptions/

```
+ /amiciOptions/
 \
  - (Attributes specifying AMICI settings, see deps/AMICI/src/hdf5.cpp)
  - sens_ind
        List of model parameter indices w.r.t. which to compute sensitivities
  - ts
        Unique list of time-points for which there is data. Determines the 
        default timepoints set to the AMICI model. 
```

## /fixedParameters/

```
+ /fixedParameters/
 \
  - conditionNames
        Any labels for columns in `/fixedParameters/k`.
        [string: nConditionVectors]
  - k
        [double: nk * nConditionVectors]
        Fixed parameter vectors for simulation or preequilibration. 
        See `/fixedParameters/simulationConditions`
  - parameterNames
        IDs of AMICI model fixed parameters
  - simulationConditions
        List of (fixed_par_preequilibration, fixed_par_simulation) tuples,
        referring to columns in `k`.
        Number of rows determines number of model simulations per cost function
        evaluation. 
        For no preequilibration, the respective values should be -1.
```

## /measurements/

```
+ /measurements/
 \
  - observableNames
        [string: ny]
        observable IDs as provided by the AMICI model
  
  Measurement-related data, each group contains datasets named 
  1..nSimulationConditions
  - t/
    [double n_t]
    Vector with timepoints for the respective rows in `y/` and `ysigma/`.
    Ascending order, maybe non-unique to specify replicates.
  - y/
  - ysigma/
    [double n_t x n_y]
    Measured values and standard deviation for each observable and timepoint.
    If sigma is given by a model parameter, entry must be NaN
```

## /parameters/

```
+ /parameters/
 \
  - lowerBound
  - upperBound
        [double n_opt_par]
        Bounds for optimization parameters. Same length and order as 
        modelParameterNames. Given in scale `pscale`. 
  - modelParameterNames
        [string n_opt_par]
        AMICI model parameter IDs
  - nominalValues
        [double n_opt_par]
        Some specified parameter vector for testing
  - optimizationSimulationMapping
        [int np * n_simulation_conditions]
        Mapping of optimization parameters to model parameters for a given 
        condition. -1 for no mapping, in which case the respective value 
        from `parameterOverrides` will be used.
  - parameterOverrides
        [int np * n_simulation_conditions]
        Constant condition-specific parameter overrides 
        (see `optimizationSimulationMapping`)
  - parameterNames
        [string n_opt_par]
        Optimization parameter names
  - pscaleOptimization
        [int n_opt_par]
        Parameter scaling for optimization parameters
  - pscaleSimulation
        [int n_simulation_conditions * n_sim_par]
        Parameter scale for each of `np` model parameters for the 
        respective condition
        (amici::ParameterScaling)

```

## Optimization options

```
+ /optimizationOptions/
 \
  - Attributes:
    - maxIter [int]
    - hierarchicalOptimization [int]
    - numStarts [int]
    - optimizer [int]
  - randomStarts [double np_opt x n_starts]

  Groups with attributes for optimizer-specific settings:
  - ceres/
  - fmincon/
  - ipopt/
  - toms611/
  See the respective source files in parpeoptimization/src/ for details.
```

## Hierarchical optimization

Those datasets are only required for use with hierarchical optimization.
If any dataset would be empty, it can be omitted.

```
/offsetParametersMapToObservables
/offsetParameterIndices
/scalingParametersMapToObservables
/scalingParameterIndices
/sigmaParametersMapToObservables
/sigmaParameterIndices
```

### (offset|scaling|sigma)ParameterIndices

1D integer array of indices to `/parameters/parameterNames`, indicating which
parameters should be considered for hierarchical optimization.
In ascending order.

### (offset|scaling|sigma)ParametersMapToObservables

Mapping of measurements to analytically computed parameters.

Integer matrix with three columns:

| parameterIdx | conditionIdx | observableIdx | 
|---|------|---|
|...|...|...|

`parameterIdx` is an index to `/(offset|scaling|sigma)ParameterIndices`

`conditionIdx` is an index to a row in `/fixedParameters/simulationConditions`

`observableIdx` is an index to `/measurements/observableNames`

Rows may occur in arbitrary order.

Different parameters for different timepoints are currently not supported.
A work-around would be introducing an additional observable.

## Simulation / training data

TODO

## Result files

TODO
