#!/usr/bin/env python3
"""
Merge parPE result files.

2017 Daniel Weindl <daniel.weindl@helmholtz-muenchen.de>

"""

import numpy as np
import h5py
import sys
import re
import os

def printUsage():
    print('Usage: %s INFILE ... OUTFILE' % __file__)
    
def mergeFiles(infiles, outfile):
    """
    Read final multistart optimization results from HDF5 infile and write to HDF5 outfile. 
    """
    
    fOut = h5py.File(outfile, "w")
    
    for infile in infiles:
        print(infile)
        with h5py.File(infile, "r") as fIn:
            copyRecursive(fIn["/"], fOut["/"])
            
def copyRecursive(inObj, rootOut):
    """
    Recursively copy inObj to the same path under rootOut
    
    NOTE: no attributes are copied for now
    """
    print("\tCopy %s to %s" % (inObj.name, rootOut.name))
    dest = rootOut.name + inObj.parent.name
    if(not inObj.name in rootOut):
        print("\tDestination %s does not exists: copy %s " %(dest, inObj.name))
        inObj.copy(inObj.name, rootOut.file[dest])
    else:
        print("\tDestination %s exists: copy members of %s" %(dest, inObj.name))
        for obj in inObj:
            print("\t\t" + obj)
            copyRecursive(inObj[obj], rootOut)

def main():
    numRequiredArgs = 2
    if len(sys.argv) <= numRequiredArgs + 1:
        printUsage()
        sys.exit(1)
        
    infiles = sys.argv[1:(len(sys.argv) - 1)]
    outfile = sys.argv[-1]
    
    mergeFiles(infiles, outfile)
    
if __name__ == "__main__":
    main()
