#!/usr/bin/env python3

import numpy as np
from math import nan
import unittest
import h5py


def checkMapping(mapping):
    # ensure consecutive scaling factor indices
    uni = np.unique(mapping[:, 0])
    assert(len(uni) == max(uni) + 1)
    numScalings = len(uni)

    # TODO: ensure non-overlapping scaling factors
    uni = np.unique(mapping[:, 1:3], axis=0)
    assert(uni.shape[0] == mapping.shape[0])


class TestHierarchical(unittest.TestCase):

    def test_okay(self):
        # 2 scaling factors, 1. for condition 1 and 2 observable 0,
        # 2. for condition 1,2,3, observable 1
        mapping = np.array([[0, 1, 0],
                            [0, 2, 0],
                            [1, 1, 1],
                            [1, 2, 1],
                            [1, 3, 1]])
        checkMapping(mapping)

    def test_overlapping(self):
        # There is no scaling factor #2 but #3
        mapping = np.array([[0, 1, 0],
                            [0, 2, 0],
                            [1, 1, 1],
                            [1, 2, 1],
                            [3, 1, 1]])
        with self.assertRaises(AssertionError):
            checkMapping(mapping)

    def test_noncontinuous(self):
        # Overlapping scaling factors 1 and 2
        mapping = np.array([[0, 1, 0],
                            [0, 2, 0],
                            [1, 1, 1],
                            [1, 2, 1],
                            [2, 1, 1]])
        with self.assertRaises(AssertionError):
            checkMapping(mapping)


def main():
    mappingToObservable = np.array([[0, 1, 0],
                                    [0, 2, 0],
                                    [1, 1, 1],
                                    [1, 2, 1],
                                    [1, 3, 1]])
    mappingToFullParams = [0, 1]
    #print(mappingToObservable)

    checkMapping(mappingToObservable)
    assert(len(mappingToFullParams)
           == len(np.unique(mappingToObservable[:, 0])))

    # used by hierarchicalOptimizationTest.cpp
    fileNameH5 = 'testhierarchical.h5'
    f = h5py.File(fileNameH5, "w")
    dset = f.create_dataset("/scalingParametersMapToObservables",
                            mappingToObservable.shape, dtype='<i4')
    dset[:] = mappingToObservable
    dset.attrs['numScalings'] = len(mappingToFullParams)

    dset = f.create_dataset("/scalingParameterIndices",
                            (len(mappingToFullParams),), dtype='<i4')
    dset[:] = mappingToFullParams


if __name__ == '__main__':
    main()
    unittest.main()
