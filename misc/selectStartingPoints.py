#!/usr/bin/env python3

"""
Script for modifying parPE optimization starting points in an HDF5 file

Creates a backup of the starting point dataset, adjusts starting point order and the number of starts. Switches off "retry" options.

2018 Daniel Weindl
"""

import h5py
import sys
import re
from typing import List

startingPointPath = '/optimizationOptions/randomStarts'
startingPointBackupPath = '/optimizationOptions/randomStartsBackup'

def printStarts(filename):
    """
    Print info on existing starting points
    """

    with h5py.File(filename, "r") as f:
        try:
            starts = f[startingPointPath]
            numStartsSelected = f['/optimizationOptions/'].attrs['numStarts']
            print("Currently selected %d starts from provided %d starting points for %d parameters" % (numStartsSelected, starts.shape[1], starts.shape[0]))

        except KeyError:
            print('{filename} does not contain starting points in {startingPointPath}.'.format(filename=filename,
                                                                                                     startingPointPath=startingPointPath))
        if startingPointBackupPath in f:
            starts = f[startingPointBackupPath]
            print("File contains starting point backup with %d starting points for %d parameters" % (starts.shape[1], starts.shape[0]))

def setStarts(filename, selection):
    """
    Create new set of starting points based on the user-selection and original starting point dataset
    """

    starts = parse_selection(selection)
    with h5py.File(filename, "r+") as f:
        backupStartsIfNotExists(f)
        print('Selecting starting points', starts)
        if startingPointPath in f: del f[startingPointPath]
        backupDset = f[startingPointBackupPath]
        f.require_dataset(startingPointPath, shape=(backupDset.shape[0], len(starts)), dtype='f8', data=backupDset[:,starts])
        updateOptions(f, len(starts))


def parse_selection(selection_str: str) -> List[int]:
    """
    Parse comma-separated list of integer ranges, return selected indices as
    integer list

    Valid input e.g.: 1 1,3 -3,4,6-7
    """
    indices = []
    for group in selection_str.split(','):
        if not re.match(r'^(?:-?\d+)|(?:\d+(?:-\d+))$', group):
            print("Invalid selection", group)
            sys.exit()
        spl = group.split('-')
        if len(spl) == 1:
            indices.append(int(spl[0]))
        elif len(spl) == 2:
            begin = int(spl[0]) if spl[0] else 0
            end = int(spl[1])
            indices.extend(range(begin, end + 1))
    return indices


def backupStartsIfNotExists(f):
    if not startingPointBackupPath in f:
        print("Creating starting point backup in", startingPointBackupPath)
        f.create_dataset(startingPointBackupPath, data=f[startingPointPath])
    else:
        print("Starting point backup exists in", startingPointBackupPath)


def updateOptions(f, numStarts):
    o = getOptionsObject(f)
    if 'retryOptimization' in o.attrs and o.attrs['retryOptimization'] != 0:
        print("Disabling 'retryOptimization'")
        o.attrs['retryOptimization'] = 0
    o.attrs['numStarts'] = numStarts


def getOptionsObject(f):
    if "/optimizationOptions" in f:
        options = f["/optimizationOptions"]
        return options

    print("Error: file does not contain /optimizationOptions")

    return None


def printUsage():
    print("""Usage: 
    selectStartingPoints.py                      -> print usage
    selectStartingPoints.py $hdf5file            -> print file info
    selectStartingPoints.py $hdf5file $selection -> sets the starting points to the selected indices
                                                    $selection: e.g.: 0,3,4-6
    """)


if __name__ == "__main__":

    if len(sys.argv) < 2 or len(sys.argv) > 3:
        printUsage()
        exit()

    filename = sys.argv[1]

    if len(sys.argv) == 2:
        printStarts(filename)
    elif len(sys.argv) == 3:
        selection = sys.argv[2]
        setStarts(filename, selection)
