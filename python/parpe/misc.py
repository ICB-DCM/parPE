import h5py
import numpy as np

"""
TODO test    
    import tempfile
    tempfile.mktemp()
"""

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
            ms = 0
            for mspath in f['/multistarts/']:
                trajectory = np.transpose(f['/multistarts/%s/iterCostFunCost'%mspath][:])
                trajectories = concatenateStarts(trajectories, trajectory)
                print(trajectory.shape, trajectories.shape)

        return trajectories


def concatenateStarts(a, b):
    """Concatenate two optimization trajectory matrices (numIteration x numStarts).
    Dimensions can be different for a and b.
    """
    if not a:
        return b

    if not b:
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
