#!/usr/bin/env python3

"""
Script for modifying parPE optimization options in an HDF5 file

Daniel Weindl 2017
"""

import h5py
import sys

def getOptionsDataset(f):
    if "/optimizationOptions" in f:
        options = f["/optimizationOptions"]
        return options
    
    print("Error: file does not contain /optimizationOptions")
    
    return None


def setOption(filename, option, value):
    f = h5py.File(filename, "r+")
    options = getOptionsDataset(f)

    if options:
        options.attrs[option] = int(value)
    else:
        exit(1)


def printOptions(filename):
    f = h5py.File(filename, "r")
    
    options = getOptionsDataset(f)
    if options:
        for o in options.attrs:
            print("%20s %3s" % (o, options.attrs[o]))
    else:
        exit(1)    


def unsetOption(filename, option):
    f = h5py.File(filename, "r+")
    
    options = getOptionsDataset(f)
    if options:
        try:
            options.attrs.__delitem__(option)
        except KeyError:
            pass
    else:
        exit(1)    


def printUsage():
    print("""Usage: 
    optimizationOptions.py              -> print All
    optimizationOptions.py key          -> print key
    optimizationOptions.py -s key value -> set key to value
    optimizationOptions.py -s key       -> remove key
    
    Options currently supported are 
              numStarts
                maxIter
              optimizer
       retryOptimization""")


if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        printUsage()
        exit()

    filename = sys.argv[1]

    if len(sys.argv) == 2:
        printOptions(filename)
        exit()
    
    if sys.argv[2] != "-s":
        print("Unrecognized options " + sys.argv(2))
        exit(1)
    
    if len(sys.argv) == 4:
        unsetOption(filename, sys.argv[3])
    
    elif len(sys.argv) == 5:
        setOption(filename, sys.argv[3], sys.argv[4])
        
