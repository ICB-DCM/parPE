#!/usr/bin/env python3
"""
Parse AMICI syms file
2017 Daniel Weindl <daniel.weindl@helmholtz-muenchen.de>

"""

import numpy as np
import pandas as pd
from libsbml import SBMLReader
import h5py
import sys
import re
import os

class AmiciSyms:
    symsFileName = ''
    
    def __init__(self, symsFileName):
        self.symsFileName = symsFileName
    
    def readFixedParameterNames(self):
        """
        Read parameter names from _syms file model.sym.k.
        """
        with open(self.symsFileName, 'r') as f:
            for l in f:
                if not l.startswith("model.sym.k"):
                    continue

                # remove "model.sym.k = [" and "]:\n"
                l = l[15:-3]
                return [s.strip() for s in l.split(",")]

    def readParameterNames(self):
        """
        Read parameter names from _syms file model.sym.p.
        """
        with open(self.symsFileName, 'r') as f:
            for l in f:
                if not l.startswith("model.sym.p"):
                    continue

                # remove "model.sym.p = [" and "]:\n"
                l = l[15:-3]
                return [s.strip() for s in l.split(",")]


        
    def readObservables(self):
        """
        Read IDs of observables as defined in the provided AMICI syms file
        """
        y = ""
        with open(self.symsFileName, 'r') as f:
            it = iter(f)
            for l in it:
                if not l.startswith("model.sym.y"):
                    continue
                y += l.strip()
                if l.find("];") >= 0:
                     break

                for l in it:
                    y += l.strip()
                    if l.find("];") >= 0:
                        break
        # remove "model.sym.y = [" and "];"
        y = y[15:-2]
        y = y.replace("...", "")
        return [o.strip() for o in y.split(",")]
    
    def readStateNames(self):
        x = ""
        with open(self.symsFileName, 'r') as f:
            it = iter(f)
            for l in it:
                if not l.startswith("model.sym.x "):
                    continue
                x += l.strip()
                if l.find("];") >= 0:
                    break
                for l in it:
                    x += l.strip()
                    if l.find("];") >= 0:
                        break
        # remove "model.sym.x = [" and "];"
        x = x[15:-2]
        x = x.replace("...", "")
        return [o.strip() for o in x.split(",")]
        

def unique(seq): 
    """
    Make unique, preserving order
    """
    seen = {}
    uni = []
    for item in seq:
        if item in seen: 
           continue
        seen[item] = True
        uni.append(item)
    return uni


def main():
    pass

if __name__ == "__main__":
    main()
    #a = AmiciSyms()
    #a.symsFileName = '/home/dweindl/src/CanPathPro-WP6-Leonard/Mouse_Models/Speedy_v3_r403445/Speedy_v3_r403445_v1_syms.m'
    #[ print(x) for x in a.readStateNames() ]
