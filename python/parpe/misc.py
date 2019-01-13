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


def getCorrTable(simulationResults, minimum_number_datapoints=0):
    """Generate correlation table by observable

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
