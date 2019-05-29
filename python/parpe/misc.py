import h5py
import numpy as np

"""
TODO test    
    import tempfile
    tempfile.mktemp()
"""

from .plotting import  correlation_coefficient


def getCostTrajectories(filename, starts = None):
    """
    Read cost trajectory from HDF5 result file.

    Arguments:
    filename: result file name
    starts: list with indices or None meaning all starts
    """
    with h5py.File(filename, 'r') as f:
        trajectories = None
        if starts:
            for ms in starts:
                trajectory = np.transpose(f['/multistarts/%s/iterCostFunCost'%ms][:])
                trajectories = concatenateStarts(trajectories, trajectory)
        else:
            for mspath in f['/multistarts/']:
                trajectory = np.transpose(f['/multistarts/%s/iterCostFunCost'%mspath][:])
                trajectories = concatenateStarts(trajectories, trajectory)

        return trajectories


def getCostTrajectory2(filename):
    """Read cost trajectory from file

    #Arguments:
    #filename: HDF5 file as generated by misc/extractMultiStartParameters.py

    Returns:
    ndarray(numIterations x numStarts): cost over iterations for all optimizations
    """

    with h5py.File(filename, 'r') as f:
        numStarts = len(f['/multistarts'])  # should take max of present ones
        numIter = 0
        # find max iter
        for ms in range(numStarts):
            ms = int(ms)
            numIter = np.max(
                [len(f['/multistarts/%d/iteration' % ms]), numIter])
        costTrajectory = np.full(shape=(numIter, numStarts), fill_value=np.nan)

        for ms in range(numStarts):
            for iteration in range(numIter):
                try:
                    if not '/multistarts/%d/iteration/%d/costFunCost' % (
                    ms, iteration) in f:
                        continue
                    value = np.nanmin(f[
                                          '/multistarts/%d/iteration/%d/costFunCost' % (
                                          ms, iteration)])
                    costTrajectory[iteration, ms] = value
                except Exception as e:
                    print(e)

    return costTrajectory


def getCostTrajectoryFromSummary(filename):
    """Read cost trajectory from file

    Arguments:
    filename: HDF5 file as generated by misc/extractMultiStartParameters.py

    Returns:
    ndarray(numIterations x numStarts): cost over iterations for all optimizations
    """

    with h5py.File(filename, 'r') as pe_summary:
        costTrajectory = pe_summary['/costTrajectory'][:]

    return costTrajectory


def concatenateStarts(a, b):
    """Concatenate two optimization trajectory matrices (numIteration x numStarts).
    Dimensions can be different for a and b.
    """
    if a is None or not a.size:
        return b

    if b is None or not b.size:
        return a

    maxIter = max(a.shape[0], b.shape[0])
    a = np.pad(a, pad_width=[(0, maxIter - a.shape[0]), (0, 0)], mode='constant', constant_values=np.nan)
    b = np.pad(b, pad_width=[(0, maxIter - b.shape[0]), (0, 0)], mode='constant', constant_values=np.nan)
    concat = np.concatenate((a, b), axis=1)

    return concat


def readSimulationsFromFile(filename):
    """
    Read result from simulations at final point from data file for all starts

    Returns:
        (measured, simulated, time, llh)[startString]
            [nCondition][nTimepoints, nObservables]
    """

    sim = {}
    mes = {}
    llh = {}
    time = {}

    with h5py.File(filename, 'r') as f:
        for ms in f['/multistarts']:
            llh[ms] = f[f'/multistarts/{ms}/llh'][:]
            mes[ms] = []
            sim[ms] = []
            time[ms] = []
            for condition_idx in range(len(f[f'/multistarts/{ms}/t'])):
                mes[ms].append(f[f'/multistarts/{ms}/yMes/{condition_idx}'][:])
                sim[ms].append(f[f'/multistarts/{ms}/ySim/{condition_idx}'][:])
                time[ms].append(f[f'/multistarts/{ms}/t/{condition_idx}'][:])

    return mes, sim, time, llh


def readSimulationsFromFileLegacy(filename):
    """
    Read result from simulations at final point from data file for all starts
    Returns:
    (measure, simulated)[startString][nCondition, nTimepoints, nObservables]
    """

    sim = {}
    mes = {}
    llh = {}

    with h5py.File(filename, 'r') as f:
        for ms in f['/multistarts']:
            mes[ms] = f['/multistarts/%s/yMes'%ms][:]
            sim[ms] = f['/multistarts/%s/ySim'%ms][:]
            llh[ms] = f['/multistarts/%s/llh'%ms][:]

    return mes, sim, llh


def getConditionNames(filename):
    with h5py.File(filename, 'r') as f:
        if '/inputData/fixedParameters/conditionNames' in f:
            conditionNames = f['/inputData/fixedParameters/conditionNames'][:]
        else:
            conditionNames = f['/fixedParameters/conditionNames'][:]

    return conditionNames


def getFinalCostFromSummary(filename):
    """Read final cost from file

    Arguments:
    filename: HDF5 file as generated by misc/extractMultiStartParameters.py

    Returns:
    ndarray(numStarts): cost at last iteration for all optimizations
    """

    with h5py.File(filename, 'r') as pe_summary:
        finalCost = pe_summary['/finalCost'][:]

    return finalCost


def getCorrTable(mes, sim, minimum_number_datapoints=0):
    """Generate correlation table by observable

    Arguments:
        simulationResults: simulationResults as obtained from
        parpe.readSimulationsFromFile

    Returns:
        ndarray with correlations of measurement and simulation for
        each start and observables in simulationResults
    """

    numStarts = len(mes)
    numCond = len(list(mes.values())[0])
    numObs = list(mes.values())[0][0].shape[1]

    corr = np.full((numStarts, numObs), np.nan)
    for ims, ms in enumerate(mes):
        ymes = None
        ysim = None

        for icond in range(numCond):
            if ymes is None:
                ymes = mes[ms][icond]
                ysim = sim[ms][icond]
            else:
                ymes = np.append(ymes, mes[ms][icond], axis=0)
                ysim = np.append(ysim, sim[ms][icond], axis=0)

        for iobservable in range(numObs):
            if np.sum(np.isfinite(
                    ymes[:, iobservable])) < minimum_number_datapoints:
                corr[ims, iobservable] = np.nan
                continue

            corr[ims, iobservable] = correlation_coefficient(
                ymes[:, iobservable], ysim[:, iobservable])

    return corr


def getCorrTableLegacy(simulationResults, minimum_number_datapoints=0):
    """Generate correlation table by observable

    Old data format

    Arguments:
        simulationResults: simulationResults as obtained from
        parpe.readSimulationsFromFile

    Returns:
        ndarray with correlations of measurement and simulation for
        each start and observables in simulationResults
    """
    numStarts = len(simulationResults[0])
    numObs = list(simulationResults[0].values())[0].shape[2]
    corr = np.full((numStarts, numObs), np.nan)
    for ims, ms in enumerate(simulationResults[0]):
        ymes = simulationResults[0][ms]
        ysim = simulationResults[1][ms]

        for observable in range(ymes.shape[2]):
            if np.sum(np.isfinite(
                    ymes[:, :, observable])) < minimum_number_datapoints:
                corr[ims, observable] = np.nan
                continue

            corr[ims, observable] = correlation_coefficient(
                ymes[:, :, observable], ysim[:, :, observable])

    return corr


class ParameterEstimationSummaryFile:
    """Interface to a parPE parameter estimation summary file as created
    with misc/extractMultiStartParameters.py
    """
    pass


class ParameterEstimationResultFile:
    """Interface to a parPE parameter estimation result file"""

    pass


def simulation_to_df(mes_df, sim, result_file, start, label, observable_ids):
    """Add simulation results as new column to PEtab measurement df"""

    with h5py.File(result_file, 'r') as f:
        condition_names = f['/inputData/fixedParameters/conditionNames'][:]
        simulation_conditions = f[
                                    '/inputData/fixedParameters/simulationConditions'][
                                :]
    cond_id_to_idx = {id_: idx for idx, id_ in enumerate(condition_names)}
    cond_comb_to_idx = {
    (condition_names[preeq_idx] if not np.isnan(preeq_idx) and preeq_idx >= 0 else -1.0,
     condition_names[sim_idx] if sim_idx >= 0 else -1.0): comb_idx
    for comb_idx, (preeq_idx, sim_idx) in enumerate(simulation_conditions)}
    obs_id_to_idx = {id_: idx for idx, id_ in enumerate(observable_ids)}

    mes_df['simulation_' + label] = np.nan
    time_tracker = {}
    for _, row in mes_df.iterrows():
        if np.isnan(row.preequilibrationConditionId):
            row.preequilibrationConditionId = -1
        condition_idx = cond_comb_to_idx[
            row.preequilibrationConditionId, row.simulationConditionId]
        observable_idx = obs_id_to_idx['observable_' + row.observableId]
        """ Replicate simulations will be identical, but lets consider them"""
        try:
            time_tracker[(condition_idx, observable_idx, row.time)] += 1
        except KeyError:
            time_tracker[(condition_idx, observable_idx, row.time)] = 0

        time_idx = time_tracker[(condition_idx, observable_idx, row.time)]

        # time_idx = 0
        mes_df.loc[row.name, 'simulation_' + label] = \
        sim[start][condition_idx][time_idx, observable_idx]
    return mes_df


def compare_optimization_results_to_true_parameters(filename: str):
    """Compare parameter estimates to true parameters. Print as table.

    Used in example notebooks.

    Arguments:
        filename: Parameter estimation result file name
    """
    with h5py.File(filename, 'r') as f:
        pscale = f['/inputData/parameters/pscaleOptimization'][:]
        names = f['/inputData/parameters/parameterNames'][:]
        true_parameters = f['/inputData/parameters/true_parameters'][:]
        expectedNllh = -f['/inputData/parameters/true_llh'][:]
        final_parameters = f['/multistarts/0/finalParameters'][:]
        exit_status = f['/multistarts/0/exitStatus'][:]
        final_cost = f['/multistarts/0/finalCost'][:]

    for i, p in enumerate(pscale):
        if p == 2:
            final_parameters[i] = np.power(10, final_parameters[i])

    print("#  __Exp____ __Act______ __Err______ __RelErr___ __ID_______")
    for i in range(len(true_parameters)):
        error = final_parameters[i]-true_parameters[i]
        rel_error = error / true_parameters[i]
        print('%d: %9.5f %11.5f %11.5f %11.5f %s' % (i, true_parameters[i],
                                       final_parameters[i],
                                       error, rel_error, names[i]))
    print()
    print('Status: %d' % exit_status)
    print('Cost: %f (expected: %f)' % (final_cost, expectedNllh)) # FIXME: expectedNllh not correctly written to file
