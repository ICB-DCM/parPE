#!/usr/bin/env python3
"""
Generate HDF5 file for parPE with fixed parameters and measurements for an AMICI-imported SBML model based on tables with fixed parameters and training data

2018 Daniel Weindl <daniel.weindl@helmholtz-muenchen.de>, Leonard Schmiester <leonard.schmiester@helmholtz-muenchen.de>

Usage: __file__ hdf5File sbmlModelFile symsModelFile cost_func exp_table

hdf5File:      Name of the generated HDF5 file
sbmlModelFile: The SBML model for which the data is generated (not used anymore, possibly in future to access annotations) 
symsModelFile: The AMICI syms file
cost_func:     The cost function file (PyBios) listing the experimental conditions to use
exp_table:     Data table for fixed parameters for all conditions listed in cost_func 
               (tsv: fixedParameters x conditions, row names according to model parameters, column names according to cost function) 
"""

import amiciHelper
import numpy as np
import pandas as pd
from libsbml import SBMLReader
import h5py
import sys
import re
import os
from termcolor import colored
import collections

SCALING_LOG10 = 2
SCALING_LIN = 0

class AmiciPyModelParser:
    def __init__(self, modelFolder):
        self.modelDir = modelFolder
    
    @staticmethod
    def parseHeader(filename):
        symbols= []
        with open(filename, 'r') as f:
            for line in f:
                splitted = line.split(' ')
                if len(splitted) == 3:
                    symbols.append(splitted[1])
                else:
                    print('Unexpected line format in %s: %s' %(filename, line))
        return symbols
    
    def readParameterNames(self):
        return self.parseHeader(os.path.join(self.modelDir, 'parameter.h'))
    
    def readFixedParameterNames(self):
        return self.parseHeader(os.path.join(self.modelDir, 'fixed_parameter.h'))
    
    def readObservables(self):
        return self.parseHeader(os.path.join(self.modelDir, 'observable.h'))
      
    
class HDF5DataGenerator:
    """    
    Generate HDF5 file with fixed parameters and measurements for an AMICI-imported SBML model
    """

    def __init__(self, fileNameSBML, fileNameSyms, fileMeasurements, fileFixedParameters):
        """
        fileNameSBML: filename of model SBML file
        fileNameSyms: filename of AMICI model syms file or the model folder of a python-amici generated model
        fileFixedParameters: filename with AMICI fixed parameter vectors for all conditions refered to in the measurement file 
        fileMeasurements: filename of measurement file (TODO: specify format)
        """

        # read symbols from compiled model to ensure correct ordering (might change during sbml import)
        if fileNameSyms.endswith('.m'):
            # AMICI model file interface
            self.amiciSyms = amiciHelper.AmiciSyms(fileNameSyms)
        else:
            self.amiciSyms = AmiciPyModelParser(fileNameSyms)

        # load input data 
        self.parseMeasurementFile(fileMeasurements)
        self.parseFixedParametersFile(fileFixedParameters)

        self.loadSBMLModel(fileNameSBML)

        # scriptDir = os.path.dirname(os.path.realpath(__file__))
        # hdf5 dataset compression
        self.compression = "gzip"
               

    def parseMeasurementFile(self, filename):
        """
        Read cost_fun file and determine number of conditions and timepoints
        """
        self.measurementDf = pd.read_csv(filename, 
                                  delimiter="\t", 
                                  index_col=False)

        print("Cost shape", self.measurementDf.shape)

        # Get list of used conditions
        # Must also consider conditionRef values, since pre-equilibration conditions might come without any datapoints,
        # and thus, not show up in condition-column
        self.conditions = amiciHelper.unique(self.measurementDf.loc[:, 'condition'].append(self.measurementDf.loc[:, 'conditionRef']))
        if np.nan in self.conditions:
            self.conditions.remove(np.nan) # might have been introduced with empty conditionRef fields
        self.numConditions = len(self.conditions)

        # when using adjoint sensitivities, we cannot keep inf -> constider late timepoint as steady-state
        print("Changing t = Inf to t = 1e8.")
        self.measurementDf.loc[self.measurementDf['time'] == np.inf, 'time'] = 1e8
        
        # list of timepoints
        self.timepoints = amiciHelper.unique(self.measurementDf.loc[:, 'time'])
        self.timepoints.sort()
        self.numTimepoints = len(self.timepoints)      
      
        print("Num conditions: ", self.numConditions)
        print("Num timepoints: ", self.numTimepoints, self.timepoints)


    def parseFixedParametersFile(self, filename):
        """
        Load file and select relevant conditions
        """
        
        self.fixedParametersDf = pd.read_csv(filename, delimiter="\t")
        self.fixedParametersDf.set_index(self.fixedParametersDf.columns.values[0], inplace=True)
        
        print("Fixed parameters orginal: ", self.fixedParametersDf.shape)
                
        # drop conditions that do not have measurements
        dropCols = [label for label in self.fixedParametersDf if not label in self.conditions]
        self.fixedParametersDf.drop(dropCols, 1, inplace=True)

        print("Fixed parameters usable: ", self.fixedParametersDf.shape)
        
        self.ensureObservedConditionsHaveFixedParameters()

    def ensureObservedConditionsHaveFixedParameters(self):
        """
        Ensure all conditions mentioned in measurement file are are present in fixed parameter file
        """
        cols = set(self.fixedParametersDf)
        for condition in self.conditions:
            if not condition in cols:
                print('Error: condition %s not in exp_table.' % condition)


    def loadSBMLModel(self, filename):
        """
        Load SBML model
        """
        # Note: must keep reader and document, otherwise model becomes invalid
        self.reader = SBMLReader()
        self.document = self.reader.readSBMLFromFile(filename)
        self.sbmlModel = self.document.getModel()


    def generateFile(self, fileNameH5):
        """
        Create the output file
        
        fileNameH5:   filename of HDF5 file that is to be generated
        """
        self.f = h5py.File(fileNameH5, "w")

        print("Generate parameter list...")
        self.generateParameterList()

        self.generateReferenceMap()

        print("Generating fixed parameters matrix...")
        self.generateFixedParameterMatrix()

        print("Generating measurement matrix...")
        self.generateMeasurementMatrix()

        print("Handling scaling parameters...")
        self.generateScalingParameterData()

        print("Copying default AMICI options...")
        self.copyAmiciOptions()

        print("Writing default optimization options...")
        self.writeOptimizationOptions()
        
        
    def generateParameterList(self):
        """ 
        Optimization to simulation parameter mapping. Write parameter names. 
        """
        
        # simulation parameters from model
        simulationParameterNames = self.amiciSyms.readParameterNames()
        self.writeStringArray("/parameters/modelParameterNames", simulationParameterNames)
        print("Number of simulation parameters: %d" % len(simulationParameterNames))
        
        # simulationParameterMap is OrderedDict with simulation parameter names as keys, and dicts with 
        # condition names as key and condition-specific instances as values. 
        # no instance name means name is equal to model parameter name  
        self.simulationParameterMap = collections.OrderedDict()
        for p in simulationParameterNames:
            self.simulationParameterMap[p] = {}
        
        # generate list of optimization parameters (>= simulation parameters because of condition-specific parameters)
        optimizationParameterNames = self.replaceScalingParameterNames(simulationParameterNames[:])
        print("Number of optimization parameters: %d" % len(optimizationParameterNames))
        self.writeStringArray("/parameters/parameterNames", optimizationParameterNames)
        assert(len(simulationParameterNames) <= len(optimizationParameterNames))
        self.generateSimulationToOptimizationParameterMapping(simulationParameterNames, optimizationParameterNames)        
        
        self.f.flush()
    
    
    def replaceScalingParameterNames(self, parameterNames):
        """
        Remove the "generic" names of the scaling parameter in the model and replace with the instances for the respective experiments 
        """
        
        conditionSpecificScalingsUsed = self.getUsedScalingParameters()
        
        #print(scalingsUsed)
        for currentConditionSpecificName in conditionSpecificScalingsUsed:
            currentGenericName = self.getGenericScalingParameterName(currentConditionSpecificName)
            #print(currentGenericName)
            if currentGenericName in parameterNames:
                parameterNames.remove(currentGenericName)
            
        parameterNames += conditionSpecificScalingsUsed
        return parameterNames

    def getUsedScalingParameters(self):
        """
        Get unique list of scaling parameters mentioned in measurement file
        with optimization parameter names (*not* as in sbml/amici model) 
        """
        
        # List of condition-specifc parameter names
        # NOTE: not using set() here which would scramble parameter order.
        # This allows use to keep starting points from /randomStarts
        # the same after regenerating the data file.
        conditionSpecificScalingParameterNames = []
        
        # Track if there have been condition-spefic names provided for this parameters. They need to be always or never condition-specific. 
        observableHasConditionSpecificParameters = {}
        for index, row in self.measurementDf.iterrows():
            conditionName = row['condition']
            conditionSpecificNames = self.splitScalingParameterNames(row['scalingParameter'])
            if len(conditionSpecificNames):
                # ensure this observable has not been encountered or was marked as conditonSpecific before
                if row['observable'] in observableHasConditionSpecificParameters and observableHasConditionSpecificParameters[row['observable']] == False:
                    raise ValueError("%s had condition-specific parameters before, but none are given for %s" % (row['observable'], row))
                else:
                    observableHasConditionSpecificParameters[row['observable']] = True
            else:
                if row['observable'] in observableHasConditionSpecificParameters and observableHasConditionSpecificParameters[row['observable']] == True:
                    raise ValueError("%s had no condition-specific parameters before, but has them in %s" % (row['observable'], row))
                else:
                    observableHasConditionSpecificParameters[row['observable']] = False
                    
            for currentConditionSpecificName in conditionSpecificNames:
                conditionSpecificScalingParameterNames.append(currentConditionSpecificName)
                currentGenericName = self.getGenericScalingParameterName(currentConditionSpecificName)
                self.simulationParameterMap[currentGenericName][conditionName] = currentConditionSpecificName

        conditionSpecificScalingParameterNames = amiciHelper.unique(conditionSpecificScalingParameterNames)
        return conditionSpecificScalingParameterNames

    def splitScalingParameterNames(self, names):
        if not isinstance(names, float):
            # N/A will be float
            # is single value or comma-separated list
            return names.split(",")
        return []
    
    def generateSimulationToOptimizationParameterMapping(self, simulationParameterNames, optimizationParameterNames):
        """
        Create dataset n_pararameters_simulation x n_conditions with indexes of respective parameters in pararameters_optimization
        """
        numSimulationParameters = len(simulationParameterNames)
        # use in-memory matrix, don't write every entry to file directly (super slow)
        parameterMap = np.zeros(shape=(numSimulationParameters, self.numConditions),)
        
        # create inverse mapping
        optimizationParameterNameToIndex = { name: idx for idx, name in enumerate(optimizationParameterNames) }        
        for idxSimulationParameter, simulationParameterName in enumerate(simulationParameterNames):
            
            if not len(self.simulationParameterMap[simulationParameterName]):
                # not a condition-specific parameter, set same index for all conditions
                parameterMap[idxSimulationParameter, :] = optimizationParameterNameToIndex[simulationParameterName]
            else:
                # condition-specific parameter
                for conditionIdx, conditionName in enumerate(self.conditions):
                    # this must be present if there is data for the observable using this parameter (checked above),
                    # otherwise we can set any index, because amici will not use this parameter and set the gradient to 0.0. we will use 0 here.
                    # TODO: handle better?
                    if conditionName in self.simulationParameterMap[simulationParameterName]:
                        conditionSpecificParameterName = self.simulationParameterMap[simulationParameterName][conditionName]
                        parameterMap[idxSimulationParameter, conditionIdx] = optimizationParameterNameToIndex[conditionSpecificParameterName]
                    else:
                        # This condition does not use the respective scaling parameter so there is no corresponding optimization parameter
                        # Cannot set to NaN in integer matrix. Will use 0. AMICI will not use this parameter anyways and its gradient will be 0.0 
                        parameterMap[idxSimulationParameter, conditionIdx] = 0

        self.f.require_dataset('/parameters/optimizationSimulationMapping', 
                               shape=(numSimulationParameters, self.numConditions), 
                               chunks=(numSimulationParameters, 1), 
                               dtype='<i4', fillvalue=0, compression=self.compression, data=parameterMap)

    def getOptimizationParameterNameForConditionSpecificSimulationParameter(self, conditionIdx, simulationParameterName, simulationParameterNames):
        scalingsForCurrentConditionByMeasurement = self.measurementDf.loc[self.measurementDf['condition'] == self.conditions[conditionIdx], 'scalingParameter']
        for scalingsForCurrentMeasurement in scalingsForCurrentConditionByMeasurement:
            if not isinstance(scalingsForCurrentMeasurement, float):
                # N/A will be float
                # is single value or comma-separated list
                for scaling in scalingsForCurrentMeasurement.split(","):
                    if self.getGenericScalingParameterName(scaling) == simulationParameterName:
                        return scaling
        return None
            


    def getGenericScalingParameterName(self, scaling):
        """
        For the given scaling parameter name from measurement file or optimization parameter names,
        get the parameter name as specified in the sbml/amici model.
        
        (i.e. remove suffix from given scaling parameter)
        """
        
        return "_".join(scaling.split("_")[:-1])
            
    def getGenericScalingParameterNameAndCondition(self, scaling):
        """
        For the given scaling parameter name from measurement file or optimization parameter names,
        get the parameter name as specified in the sbml/amici model.
        
        (i.e. remove suffix from given scaling parameter)
        """
        
        splitted = scaling.split("_")
        genericName = "_".join(splitted[:-1])
        suffix = splitted[-1]
        return genericName, suffix
            
    def getGenericScalingParameterNames(self, scalingsUsed):
        """
        For the given scaling parameter names from measurement file or optimization parameter names,
        get unique list of the parameter names as specified in the sbml/amici model.
        """
        # NOTE: not using set() here which would scramble parameter order.
        # This allows use to keep starting points from /randomStarts
        # the same after regenerating the data file.  
        genericNames = []
        for s in scalingsUsed:
            genericName = self.getGenericScalingParameterName(s)
            genericNames.append(genericName)
        return amiciHelper.unique(genericNames)
    
    
    def generateFixedParameterMatrix(self):
        """
        Write fixed parameters dataset. 
        """
        
        k = self.amiciSyms.readFixedParameterNames()
        self.nk = len(k)
        print("Number of fixed parameters: %d" % len(k))

        # Create in-memory table, write all at once for speed
        fixedParameterMatrix = np.full(shape=(self.nk, self.numConditions), fill_value=np.nan)
        for i in range(len(k)):
            self.handleFixedParameter(i, k[i], fixedParameterMatrix)

        dset = self.createFixedParameterDatasetAndWriteAttributes(k, fixedParameterMatrix)
        
        self.f.flush()


    def createFixedParameterDatasetAndWriteAttributes(self, fixedParameters, data):
        """
        Create fixed parameters data set and annotations
        """
        g = self.f.require_group("/fixedParameters")

        self.writeStringArray("/fixedParameters/parameterNames", fixedParameters)
        self.writeStringArray("/fixedParameters/conditionNames", self.conditions)      
              
        # chunked for reading experiment-wise
        nk = len(fixedParameters)
        dset = self.f.create_dataset("/fixedParameters/k",
                                     (nk, self.numConditions),
                                     dtype='f8', chunks=(nk, 1), compression=self.compression, data=data)

        # set dimension scales 
        dset.dims.create_scale(self.f['/fixedParameters/parameterNames'], 'parameterNames')
        dset.dims.create_scale(self.f['/fixedParameters/conditionNames'], 'conditionNames')
        dset.dims[0].attach_scale(self.f['/fixedParameters/parameterNames'])
        dset.dims[1].attach_scale(self.f['/fixedParameters/conditionNames'])
        
        return dset

    def generateReferenceMap(self):
        """
        Write index map reference conditions to be used for pre-equilibration
        
        TODO self.measurementDf.conditionRef: this should go to exp table to avoid potential redundancy/ambiguity in measurement file
        """
        
        referenceMap = list(range(len(self.conditions))) # default: no other reference conditions
        for index, row in self.measurementDf.iterrows():
            conditionIdx = self.conditions.index(row['condition'])
            if isinstance(row['conditionRef'], float) and  np.isnan(row['conditionRef']):
                conditionIdxRef = conditionIdx 
            else:
                conditionIdxRef = self.conditions.index(row['conditionRef'])
            
            if referenceMap[conditionIdx] != conditionIdxRef and referenceMap[conditionIdx] != conditionIdx:
                printf("Error: Ambiguous assignment of reference conditions for %s: %s vs. %s" 
                       % (self.conditions[conditionIdx], self.conditions[referenceMap[conditionIdx]], self.conditions[conditionIdxRef])) 
            
            referenceMap[conditionIdx] = conditionIdxRef
        dset = self.f.create_dataset("/fixedParameters/referenceConditions", shape=(self.numConditions,), dtype="<i4", data=referenceMap)

    def writeStringArray(self, path, strings):
        """
        Write string array to hdf5
        """
        dt = h5py.special_dtype(vlen=str)
        dset = self.f.create_dataset(path, (len(strings),), dtype=dt)
        dset[:] = [ s.encode('utf8') for s in strings ]
        self.f.flush()


    def writeFloatArray(self, path, values):
        """
        Write float array to hdf5
        """
        dset = self.f.create_dataset(path, (len(values),), dtype='f8')
        dset[:] = values
        self.f.flush()

    def writeIntArray(self, path, values):
        """
        Write integer array to hdf5
        """
        dset = self.f.create_dataset(path, (len(values),), dtype='<i4')
        dset[:] = values
        self.f.flush()


    def handleFixedParameter(self, parameterIndex, parameterName, dset):
        """
        Extract parameter values from data table and write to HDF5 dataset.

        NOTE: parameters which are not provided in fixedParametersDf are set to 0.0
        """
        
        if parameterName in self.fixedParametersDf.index:
            #print(parameterName, ' ' , parameterIndex)
            #print(self.fixedParametersDf.loc[parameterName, self.conditions].values)
            dset[parameterIndex, :] = self.fixedParametersDf.loc[parameterName, self.conditions].values
        else:
            # initial states
            #         s = self.sbmlModel.getSpecies(parameterName)
            print("Fixed parameter not found in ExpTable, setting to 0.0: ", parameterName)
            dset[parameterIndex, :] = 0.0


    def generateMeasurementMatrix(self):
        """
        Generate matrix with training data for all observables defined in the _syms file.
        
        TODO: at the moment we cannot deal with replace measurements (same condition, same observable, same timepoint) since not supported by AMICI
        should alert user in case that happens, currently will be silently overwritten by the latter 
        dirty workaround: take mean of change time to t+eps  
        
        NOTE: data is expected to have numTimepoints as leading dimensions. parameter estimation might fail silently if not.
        """

        self.f.create_group("/measurements")
        
        self.y = self.amiciSyms.readObservables()
        self.ny = len(self.y)
        self.observables = self.getObservableNamesFromSbml()
        self.checkObservablesSbmlMatchAmici(self.observables, self.y)
        self.writeStringArray("/measurements/observableNames", self.observables)

        print("Number of observables: %d" % self.ny)

        dsetY = self.f.create_dataset("/measurements/y", 
                                      shape=(self.numConditions, self.numTimepoints, self.ny), 
                                      chunks=(1, self.numTimepoints, self.ny), 
                                      fillvalue=np.nan, dtype='f8', compression=self.compression)
        dsetSigmaY = self.f.create_dataset("/measurements/ysigma", 
                                      shape=(self.numConditions, self.numTimepoints, self.ny), 
                                           chunks=(1, self.numTimepoints, self.ny), 
                                           fillvalue=1.0, dtype='f8', compression=self.compression)     
        
        self.writeMeasurements(dsetY, dsetSigmaY)
        self.f.flush()
        
        #TODO one is measurement reference -> need for generating simulation data
        #dset = self.f.create_dataset("/fixedParameters/reference",
        #                             (nk, self.numConditions),
        #                             dtype='f8', chunks=(nk, 1))
        # TODO reference = map[condition][observable]
        # self.writeStringArray("/fixedParameters/celllineNames", self.cells)
        # self.writeIntArray("/fixedParameters/celllineIdx", self.getConditionCelllineIndexMap())
        #g.attrs['numCelllines'] = len(set(self.cells))
            

    def getObservableNamesFromSbml(self):
        """
        Get array with names of observables from SBML models.
        
        (Observables are not accounted for in the SBML standard. 
        We implement them as parameters and assigment rules.
        The respective parameters start with "observable_",
         but do not end in "_sigma", which is reserved for encoding the error model)  
        """
        observables = []
        for p in self.sbmlModel.getListOfParameters():
            p = p.getId()
            if p.startswith('observable_') and not p.endswith('_sigma'):
                #observableName = p[len('observable_'):]
                observableName = p
                observables.append(observableName)
        return observables
    
        
    def writeMeasurements(self, dsetY, dsetSigmaY):
        """
        Write measurements to hdf5 dataset
        """            
        # write data from each row in measurementDf
        for index, row in self.measurementDf.iterrows():
            conditionIdx = self.conditions.index(row['condition'])
            observableIdx = self.observables.index(row['observable'])
            timeIdx = self.timepoints.index(row['time'])
            dsetY     [conditionIdx, timeIdx, observableIdx] = float(row['measurement'])
            dsetSigmaY[conditionIdx, timeIdx, observableIdx] = float(row['sigma'])

    def checkObservablesSbmlMatchAmici(self, observables, y):
        """
        Check if order of observables in SBML model match that of AMICI model (which is assumed here)
        
        Problem: order may have changed due to SBML which does not necessarily conserve ordering
        Problem: formula matching is not straightforward, since they might have changed due to symbolic processing  
        """
        
        assert(len(y) == len(observables))

        for i in range(len(observables)):
            observableName = observables[i]
            amiciObsFormula = y[i]
            # rough check if the sbml and model observable match
            obsSpecies = observableName[(observableName.rfind('_') + 1):]
            if amiciObsFormula.find(obsSpecies) < 0:
                print(colored('Warning: Not sure if SBML observable "%s" matches AMICI formula "%s"' % (obsSpecies, amiciObsFormula), 'yellow'))                  


    def generateScalingParameterData(self):
        """
        Deal with offsets, proportionality factor and sigmas for hierarchical optimization
        
        Generate the respective mappings
        """
        
        # TODO should check at the end of the function if scalingIndices lists are non-overlapping 
        self.ensureNonOverlappingParameterForHierarchicalOptimization()
        
        self.handleProportionalityFactors() 
        self.handleOffsetParameter() # must call after handleProportionalityFactors
        self.handleSigmas()


    def ensureNonOverlappingParameterForHierarchicalOptimization(self):
        """
        Current parPE implementation of hierarchical optimization assumes that only one single parameter 
        (sigma, offset or proportionality) is to be calculated per observable. Will fail silently if there are more, 
        thus we need to check here.
        """
        
        map = {}
        for index, row in self.measurementDf.iterrows():
            scalings = self.splitScalingParameterNames(row['scalingParameter'])
            if not len(scalings):
                continue
            if len(scalings) > 1:
                print("Warning: multiple scaling parameters found for\n\t%s\n\twhich one to choose for hierarchical optimization? Assuming first" % row)
            #genericScaling = self.getGenericScalingParameterName(scalings[0])         
            if row['observable'] in map and map[row['observable']] != scaling:
                print("Warning: multiple scaling parameters found for\n\t%s\n\twhich one to choose for hierarchical optimization? previously found %s" % (row, map[row['observable']]))
            
            map[row['observable']] = scaling

    
    def handleOffsetParameter(self):
        """
        Write list of offset parameters selected for hierarchical optimization
        TODO : not yet implemented
        """
        
        observablesBlacklist = set()
        if '/scalingParametersMapToObservables' in self.f:
            set(self.f["/scalingParametersMapToObservables"][:, 2])
                                      
        # TODO: ensure that this observable does not have a proportionality factor
        parameterNames = self.f['/parameters/parameterNames'][:].tolist()
        offsetsForHierarchical = [x for x in self.getUsedScalingParameters() if x.startswith("offset_") and parameterNames.index(x) not in observablesBlacklist ]
        offsetsForHierarchicalIndices = [ parameterNames.index(x) for x in offsetsForHierarchical ]
        [ self.ensureOffsetIsOffsetWithPositiveSign(x) for x in offsetsForHierarchical ]
        print("Number of offset parameters for hierarchical optimization: %d" % len(offsetsForHierarchicalIndices))
        
        if len(offsetsForHierarchicalIndices) == 0:
            return
        
        dset = self.f.require_dataset("/offsetParameterIndices", 
                                      shape=(len(offsetsForHierarchicalIndices),), 
                                      dtype='<i4', 
                                      data=offsetsForHierarchicalIndices)

        # find usages for the selected parameters
        use = []
        for index, row in self.measurementDf.iterrows():
            currentScalings = self.splitScalingParameterNames(row['scalingParameter'])
            for s in currentScalings:
                #print(s, offsetsForHierarchical)
                if s in offsetsForHierarchical:
                    scalingIdx = offsetsForHierarchical.index(s)
                    conditionIdx = self.conditions.index(row['condition'])
                    observableIdx = self.observables.index(row['observable'])
                    use.append((scalingIdx, conditionIdx, observableIdx))
       
        dset = self.f.require_dataset("/offsetParametersMapToObservables", 
                                      shape=(len(use), 3), 
                                      dtype='<i4')
        dset[:] = use
                
    def ensureOffsetIsOffsetWithPositiveSign(self, scaling):
        """
        Current parPE implementation of hierarchical optimization assumes that offset parameters have positive sign.
        Will fail silently if this is not true, therefore check here that it is in indeed an offset parameter and 
        that it has positive sign
        
        TODO: sympy
        """
        
        print(colored("Ensure that %s is selected correctly as offset for hierarchical optimization and has positive sign (%s)." 
                      % (scaling, self.getFormulaForScalingParameter(scaling)), "yellow"))
        
        
    def getFormulaForScalingParameter(self, scaling):
        """
        Get formula for the (first) observable associated with the provided scaling factor 
        (name provided as in measurement list, not as in model) 
        """
        name = self.getGenericScalingParameterName(scaling)
        formula = ''
        for i in range(len(self.y)):
            if self.y[i].find(name) > -1:
                formula = self.y[i]
                break
        return formula


    def handleProportionalityFactors(self):
        """
        Write datasets specifying which parameters to consider for hierarchical optimization
        
        TODO: 
        Figure out how to deal with observables with multiple scaling parameters; which one is preferred? proportionality over offset?
        
        TODO: 
        how to let the user select which ones to add?
        """        
       
        parameterNames = self.f['/parameters/parameterNames'][:].tolist()

        scalingsForHierarchical = [x for x in self.getUsedScalingParameters() if x.startswith("scaling_") ]
        if not len(scalingsForHierarchical):
            return
        scalingsForHierarchicalIndices = [ parameterNames.index(x) for x in scalingsForHierarchical ]
        [ self.ensureProportionalityFactorIsProportionalityFactor(x) for x in scalingsForHierarchical ]

        dset = self.f.require_dataset("/scalingParameterIndices", 
                                      shape=(len(scalingsForHierarchicalIndices),), 
                                      dtype='<i4', 
                                      data=scalingsForHierarchicalIndices)
        print("Number of proportionality factors for hierarchical optimization: %d" % len(scalingsForHierarchicalIndices))
       
        # find usages for the selected parameters
        use = []
        for index, row in self.measurementDf.iterrows():
            currentScalings = self.splitScalingParameterNames(row['scalingParameter'])
            for s in currentScalings:
                #print(s, scalingsForHierarchical)
                if s in scalingsForHierarchical:
                    scalingIdx = scalingsForHierarchical.index(s)
                    conditionIdx = self.conditions.index(row['condition'])
                    observableIdx = self.observables.index(row['observable'])
                    use.append((scalingIdx, conditionIdx, observableIdx))
       
        dset = self.f.require_dataset("/scalingParametersMapToObservables", 
                                      shape=(len(use), 3), 
                                      dtype='<i4')
        dset[:] = use
        

    def ensureProportionalityFactorIsProportionalityFactor(self, scaling):
        """
        Ensure that this is a proportionality factor (a as in y = ax + b)
        
        TODO sympy
        """
    
        formula = self.getFormulaForScalingParameter(scaling)
        print(colored("Ensure that %s is selected correctly as proportionality factor for hierarchical optimization (%s)." % (scaling, formula), "yellow"))
    
    
    def handleSigmas(self):
        """
        Deal with sigma parameters
        
        TODO : not yet implemented
        """
        parameterNames = self.f['/parameters/parameterNames'][:].tolist()
        sigmasForHierarchical = [x for x in self.getUsedScalingParameters() if x.startswith("_sigma") ]
        
        if len(sigmasForHierarchical):
            print(colored("Sigmas currently not supported (%s)." % (sigmasForHierarchical), "yellow"))

    
    def copyAmiciOptions(self):
        """
        Write simulation options
        """
        g = self.f.require_group("/amiciOptions")
        #g.attrs['qpositivex'] = [0.0] * len(self.amiciSyms.readStateNames())
        #g.attrs['kappa'] = [np.nan] * len(self.amiciSyms.readParameterNames())
        #g.attrs['theta'] = [np.nan] * self.nk
        g.attrs['sensi'] = 1
        g.attrs['sensi_meth'] = 2
        g.attrs['tstart'] = 0.0
        g.attrs['atol'] = 1e-16
        g.attrs['interpType'] = 1
        g.attrs['ism'] = 1
        g.attrs['iter'] = 2
        g.attrs['linsol'] = 9
        g.attrs['lmm'] = 2
        g.attrs['maxsteps'] = 10000
        g.attrs['newton_preeq'] = 0
        g.attrs['nmaxevent'] = 0
        g.attrs['ordering'] = 0
        g.attrs['rtol'] = 1e-8
        g.attrs['stldet'] = 1

        np = self.f['/parameters/modelParameterNames'].shape[0]

        self.f.require_dataset('/amiciOptions/sens_ind', shape=(np,), dtype="<i4", data=range(np))

        self.f.require_dataset('/amiciOptions/ts', shape=(len(self.timepoints),), dtype="f8", data=self.timepoints)

        # set pscale based on whether is scaling parameter (log10 for non-hierarchical, lin for hierarchical)
        linParametersAmiciIndices = self.getAnalyticallyComputedSimulationParameterIndices()
        self.f.require_dataset('/amiciOptions/pscale', shape=(np,), dtype="<i4", data=[2 * (ip not in linParametersAmiciIndices) for ip in range(np) ])
        
        # TODO: move out of here 
        np = self.f['/parameters/parameterNames'].shape[0]
        linParametersOptimizationIndices = self.getAnalyticallyComputedOptimizationParameterIndices()
        self.f.require_dataset('/parameters/parameterScaling', shape=(np,), dtype="<i4", data=[2 * (ip not in linParametersOptimizationIndices) for ip in range(np) ])
        # TODO for above parameters, the bounds must be adapted for non-hierarchical optimization!
        #  or use different pscale then
        
    def getAnalyticallyComputedSimulationParameterIndices(self):
        """
        Get model parameter index (not optimization parameter index) of all analytically computed parameters
        """
        parameterNamesModel = []
        if '/offsetParameterIndices' in self.f:
            parameterNamesOptimization = self.f['/parameters/parameterNames'][self.f['/offsetParameterIndices']]
            parameterNamesModel.extend(set([self.getGenericScalingParameterName(o) for o in parameterNamesOptimization]))
            
        if '/scalingParameterIndices' in self.f:
            parameterNamesOptimization = self.f['/parameters/parameterNames'][self.f['/scalingParameterIndices']]
            parameterNamesModel.extend(set([self.getGenericScalingParameterName(o) for o in parameterNamesOptimization]))

        if '/sigmaParameterIndices' in self.f:
            parameterNamesOptimization = self.f['/parameters/parameterNames'][self.f['/sigmaParameterIndices']]
            parameterNamesModel.extend(set([self.getGenericScalingParameterName(o) for o in parameterNamesOptimization]))

        return [ self.f['/parameters/modelParameterNames'][:].tolist().index(p) for p in set(parameterNamesModel) ]
    
    
    def getAnalyticallyComputedOptimizationParameterIndices(self):
        """
        Get optimization parameter index of all analytically computed parameters
        """
        indices = []
        if '/offsetParameterIndices' in self.f:
            indices.extend(self.f['/offsetParameterIndices'])
            
        if '/scalingParameterIndices' in self.f:
            indices.extend(self.f['/scalingParameterIndices'])

        if '/sigmaParameterIndices' in self.f:
            indices.extend(self.f['/sigmaParameterIndices'])

        return list(set(indices))
        

    def writeOptimizationOptions(self):
        """
        Create groups and write some default optimization settings
        """
        
        # set common options
        g = self.f.require_group('optimizationOptions')
        g.attrs['optimizer'] = 0 # IpOpt
        g.attrs['retryOptimization'] = 1
        g.attrs['hierarchicalOptimization'] = 1
        g.attrs['numStarts'] = 1

        # set IpOpt options
        g = self.f.require_group('optimizationOptions/ipopt')
        g.attrs['max_iter'] = 100
        g.attrs['hessian_approximation'] = np.string_("limited-memory")
        g.attrs["limited_memory_update_type"] = np.string_("bfgs")
        g.attrs["tol"] = 1e-9
        g.attrs["acceptable_iter"] = 1
        g.attrs["acceptable_tol"] = 1e20 # set ridiculously high, so only the acceptable_* options below matter
        g.attrs["acceptable_obj_change_tol"] = 1e-18
        g.attrs["watchdog_shortened_iter_trigger"] = 0;
                
        # set CERES options
        g = self.f.require_group('optimizationOptions/ceres')
        g.attrs['max_num_iterations'] = 100
        
        # set toms611/SUMSL options
        g = self.f.require_group('optimizationOptions/toms611')
        g.attrs['mxfcal'] = 1e8
    
        self.writeBounds()
        self.writeStartingPoints()
        
        
    def writeBounds(self):
        """
        Parameter bounds for optimizer
        """
        numParams = self.f['/parameters/parameterNames'].shape[0]
        min = self.f.require_dataset('/parameters/lowerBound', [numParams], 'f8')
        max = self.f.require_dataset('/parameters/upperBound', [numParams], 'f8')
        min[:] = [-2] * numParams
        max[:] = [2] * numParams


    def writeStartingPoints(self):
        """
        Write a list of random starting points uniformly sampled from the parameter bounds. 
        Parameter bounds need to be written beforehand.
        """  
        numParams = self.f['/parameters/parameterNames'].shape[0]
        numStartingPoints = 100
        np.random.seed(0)
        startingPoints = self.f.require_dataset('/optimizationOptions/randomStarts', [numParams, numStartingPoints], 'f8')
        min = self.f['/parameters/lowerBound'][:]
        max = self.f['/parameters/upperBound'][:]
        startingPoints[:] = np.transpose(np.transpose(np.random.rand(numParams, numStartingPoints)) * (max-min) + min)  


def main():  
    if len(sys.argv) < 6:
        print("Usage: %s hdf5File sbmlModelFile symsModelFile/pyModelFolder cost_func exp_table" % __file__)
        sys.exit()

    (hdf5File, sbmlModelFile, symsModelFile, costFunctionFile, expTableFile) = sys.argv[1:6]
    h5gen = HDF5DataGenerator(sbmlModelFile, symsModelFile, costFunctionFile, expTableFile)
    h5gen.generateFile(hdf5File)

if __name__ == "__main__":
    main()
