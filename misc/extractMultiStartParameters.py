#!/usr/bin/env python3
"""
Extract multi start results from parameter estimation result files.

2017 Daniel Weindl <daniel.weindl@helmholtz-muenchen.de>

"""

import numpy as np
import h5py
import sys
import re
import os

def printUsage():
    print('Usage: %s INFILE OUTFILE' % __file__)
    
def extractMultiStartResults(infile, outfile):
    """
    Read final multistart optimization results from HDF5 infile and write to HDF5 outfile. 
    """
    
    fIn = h5py.File(infile, "r")
    numParameters = len(fIn['/multistarts/1/iteration/0/costFunParameters'])
        
    msGroup = fIn['/multistarts']
    numStarts = len(msGroup)

    fOut = h5py.File(outfile, "w")
    
    dsetParameters = fOut.create_dataset("/finalParameters",
                                 (numParameters, numStarts),
                                  dtype='f8', fillvalue=np.nan)
    
    dsetParameterTrajectories = fOut.create_dataset("/parameterTrajectory",
                                 (numParameters, 0, numStarts),
                                 maxshape=(numParameters, None, numStarts),
                                 chunks=True,
                                 dtype='f8', fillvalue=np.nan)

    dsetCostTrajectories = fOut.create_dataset("/costTrajectory",
                                 (0, numStarts),
                                 maxshape=(None, numStarts),
                                 chunks=True,
                                 dtype='f8', fillvalue=np.nan)

    dsetExitStatus= fOut.create_dataset("/exitStatus",
                                 (1, numStarts), dtype='<i4', fillvalue=-1)
    
    dsetFinalCost = fOut.create_dataset("/finalCost",
                                 (1, numStarts),
                                  dtype='f8', fillvalue=np.nan)
    
    for start in range(numStarts):
        print('Copying /multistarts/%d/' % start)
        try:
            # extend to max number of iterations if necessary
            trajectory = fIn['/multistarts/%d/iterCostFunParameters' % start]
            if trajectory.shape[1] > dsetParameterTrajectories.shape[1]:
                dsetParameterTrajectories.resize(trajectory.shape[1], 1)
            dsetParameterTrajectories[:, 0:trajectory.shape[1], start] = trajectory[:,:]

            trajectory = fIn['/multistarts/%d/iterCostFunCost' % start]
            if trajectory.shape[1] > dsetCostTrajectories.shape[0]:
                dsetCostTrajectories.resize(trajectory.shape[1], 0)
            dsetCostTrajectories[0:trajectory.shape[1], start] = trajectory[:]

            # try copy last, since they might not exist
            dsetExitStatus[:, start] = fIn['/multistarts/%d/exitStatus' % start]
            dsetParameters[:, start] = fIn['/multistarts/%d/finalParameters' % start]
            dsetFinalCost[:, start]  = fIn['/multistarts/%d/finalCost' % start]
        except KeyError:
            pass


    
def main():
    numRequiredArgs = 2
    if len(sys.argv) != numRequiredArgs + 1:
        printUsage()
        sys.exit(1)
        
    (infile, outfile) = sys.argv[1:(numRequiredArgs + 1)]
    
    extractMultiStartResults(infile, outfile)
    
if __name__ == "__main__":
    main()
