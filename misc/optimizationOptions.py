#!/usr/bin/env python3

"""
Script for modifying parPE optimization options in an HDF5 file

Daniel Weindl 2017
"""

import h5py
import sys
import re 
import numpy as np

def getOptionsObject(f):
    if "/optimizationOptions" in f:
        options = f["/optimizationOptions"]
        return options
    
    print("Error: file does not contain /optimizationOptions")
    
    return None


def setOption(filename, option, value):
    """
    Set the given option to value 
    """
    value = convertValue(value)
    
    f = h5py.File(filename, "r+")
    options = getOptionsObject(f)

    if options:
        parts = option.split('/')
        if len(parts) == 1:
            options.attrs[option] = value
        else:
            g = options.require_group(parts[0])
            g.attrs[parts[1]] = value
    else:
        exit(1)

def convertValue(value):
    """
    Guess type of data represented as string value and convert to the respective type 
    """
    
    if re.match(r'^\d+$', value):
        return int(value)
    
    try:
        return float(value)
    except ValueError:
        return np.string_(value)

def printOptions(filename):
    """
    Recursively print all options in /optimizationOptions
    """
    f = h5py.File(filename, "r")
    
    options = getOptionsObject(f)
    if options:
        printAttributes(options)
      
        for o in options:
            printAttributes(options[o], o + "/")

    else:
        exit(1)    

def printAttributes(object, prefix=''):
    """
    Print all attributes of the given object 
    """
    
    for o in object.attrs:
        print("%40s %12s" % (prefix + o, object.attrs[o]))
    
    
def unsetOption(filename, option):
    f = h5py.File(filename, "r+")
    
    group = getOptionsObject(f)
    if group:
        parts = option.split('/')
        if len(parts) > 1:
            group = group.require_group(parts[0])
            option = parts[1]
        try:
            group.attrs.__delitem__(option)
        except KeyError:
            pass
    else:
        exit(1)    


def printUsage():
    print("""Usage: 
    optimizationOptions.py              -> print All
    optimizationOptions.py key          -> print key
    optimizationOptions.py -s key value -> set key to value (use e.g. optimizername/option for optimizer-specific values)
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
        
