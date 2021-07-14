# Mini-batch optimization with parPE

Mini-batch optimization usually refers to gradient-based local optimization 
strategies such as stochastic gradient descent (SGD) with batch sizes between
1 (sometimes referred to as online learning) and the whole data set (full batch
optimization). In the context of parameter estimation for ODE models, such as 
done in parPE, this means that in each step of optimization only a subset of 
all experimental conditions is simulated. From those simulations, an estimate 
of the objective function and its gradient is computed and used for 
local optimization. ParPE has its own mini-batch optimization solver 
implemented, which includes a series of different mini-batch optimization 
algorithms, tuning possibilities and improvement tailored to parameter 
estimation of ODE models. More background on the method can be found int the
[following preprint](https://www.biorxiv.org/content/10.1101/859884v1).


## Using the mini-batch optimization solver of parPE

For the user of parPE, mini-batch optimization works similar to full-batch 
optimization with Ipopt or Fides, only a group of optimizer settings has to be
changed. For this purpose, we assume that the folder `parPE/misc` (which 
includes the file `optimizationOptions.py`) is in the python path.
1. The optimizer itself has to be changed to the mini-batch solver via \
`optimizationOptions.py <name_of_hdf5_input_file.h5> -S optimizer 10`
2. At the moment, hierarchical optimization has to be disabled, as its 
incompatible with the current way how mini-batches are chosen: \
`optimizationOptions.py <name_of_hdf5_input_file.h5> -S hierarchicalOptimization 0`
3. As the computation time (but not necessarily the wall time) is typically 
reduced when using mini-batch optimization, more local optimizations can be run
for most problems. Those can be set via 
`optimizationOptions.py <name_of_hdf5_input_file.h5> -S numStarts <number_of_starts>`


## Setting mini-batch optimization specific hyperparameters

Beyind the mentioned settings, mini-batch optimizers have a couple of specific
hyperparameters which need to be set for the mini-batch optimizer itself.
The most intuitive ones are the number of epochs, i.e., passes through the
whole dataset, which roughly corresponds to the maximum number of iterations 
in full-batch optimization, and the mini-batch size. Those can be set via \
`optimizationOptions.py <name_of_hdf5_input_file.h5> -S minibatch/batchSize 
<intended_batch_size>` 
or 
`optimizationOptions.py <name_of_hdf5_input_file.h5> -S minibatch/maxEpochs 
<number_of_epochs>`, respectively.

### Optimization algorithms
The optimizations algorithm can be set via 
`optimizationOptions.py <name_of_hdf5_input_file.h5> -S 
minibatch/parameterUpdater <name_of_algorithm>`, where `<name_of_algorithm>` 
can be chosen among the following ones:
 * "Adam": Adapted version of the Adam algorithm by Kingma et al.
 * "AdamClassic": Original version of the Adam algorithm by Kingma et al.
 * "RmsProp": RMSProp algorithm as proposed by Hinton
 * "Vanilla": Simple ("Vanilla") stochastic gradient descent
 * "Momentum": Stochastic gradient descent with (vanilla) momentum

### Learning rate scheduling
The next important hyperparameter is the learning rate, which influences the 
step size of the optimizer. ParPE allows to choose an initial learning rate,
at final learning rate (for the last epoch), and a mode of interpolation 
between those two. These things can be set via:
 * `optimizationOptions.py <name_of_hdf5_input_file.h5> -S 
minibatch/startLearningRate <floating_point_number>` 
 * `optimizationOptions.py <name_of_hdf5_input_file.h5> -S 
minibatch/endLearningRate <floating_point_number>`
 * `optimizationOptions.py <name_of_hdf5_input_file.h5> -S 
minibatch/learningRateInterpMode <interpolation_mode>`, where the interpolation
mode can be one out of `linear`, `inverseLinear`, and `logarithmic`, where 
logarithmic interpolation is the default.

### ODE-model specific implementations: Rescue interceptor and line-search
In addition to those settings, parPE implements some specific improvements for
mini-batch optimization of ODE models: First, it has a "rescueInterceptor",
which tries to recover local optimization runs if the objective function cannot
be evaluated, by shrinking the step size. It can be activated via
`optimizationOptions.py <name_of_hdf5_input_file.h5> -S 
minibatch/rescueInterceptor <integer flag>`,
where "integer flag" can be one of
 * 0 (interceptor switched off)
 * 1 (interceptor switched on)
 * 2 (interceptor switched on *and* optimization is trying a cold restart if 
 the run cannot be recovered despite reaching the maximum number of step size 
 shrinkage steps)

Second, parPE allows to perform mini-batch optimization using line-search.
Line-search can be enabled by setting the option
`optimizationOptions.py <name_of_hdf5_input_file.h5> -S 
minibatch/lineSearchSteps <maximum_number_of_steps>`, where the maximum number 
of line-search steps needs to be set. A maximum of 3 to 5 steps has shown to 
be efficient.

## Choosing hyperparameters

### The learning rate
Keep in mind, that a good learning rate depends on the optimization algorithm: 
For vanilla SGD, the learning rate is equal to the step size in parPE, whereas 
for Adam and RmsProp, the learning rate is multiplied with a more sophisticated 
parameter update, which has in general a step length larger than 1. Hence,
a good learning rate for Vanilla SGD will have larger values than for RmsProp 
or Adam (at least for the implementation in parPE).
For many applications, typical initial optimizers step sizes are often in the 
range of `[0.1, 10]`, and higher dimensional problems, step sizes are typically
larger than for lower dimensional problems.
For Adam and RmsProp, the initial parameter update, which gets multiplied with 
the square root of the problem dimension. Hence, if the problems has 10.000 
free parameters, the initial step size is a 100-fold larger than the learning 
rate. For RmsProp, the step size behaves similar to the one of Adam except that
the step sizes is (initially) biased towards 0 and hence smaller (up to at most 
1 order of magnitude). For SGD with momentum, step sizes should behave similar
to those of vanilla SGD. These things should be considered when deciding for a
learning rate schedule.

### The mini-batch size
In first tests, small mini-batch sizes worked clearly better than large ones.
At the moment, there is is no way to a priori determine the best mini-batch 
size, but this hyperparameter must be chosen in some trial-and-error manner.
Fur this purpose, it seems reasonable to start with very small mini-batch 
sizes and then gradually increase them and check how this affects the 
optimization results.

## Reporting of objective function values and gradients: Important differences 
As already mentioned, parPE only evaluates the objective function and the 
gradient on a subset of the data set in each optimization step. Hence, the 
reported values in the optimization history are only estimates based on the
current optimization step, and typically highly stochastic.
If optimization is finished successfully, the objective function is recomputed 
once on the whole data set (values stored in the last optimization step).
However, it can make sense to post-process the optimization runs, when 
assessing the optimization result: If a run has finished, the values reported 
as "finalParameter" can be used. Otherwise, it makes sense to simply use the
last reported parameter vector of the optimization history for a specific run.
These extracted final parameter can then be used to compute the actual 
objective function value on the whole data set, also including hierachically
computed nuisance parameters, if present. Such a post-processing step and 
post-optimization evaluation allows at the moment the best possible insight 
about the optimization results.
