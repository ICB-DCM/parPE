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
import math
import sympy as sp

SCALING_LOG10 = 2
SCALING_LIN = 0
# index for no reference/preequilibration condition
NO_PREEQ_CONDITION_IDX = -1


class AmiciPyModelParser:
    """Parser for AMICI-generated header files"""

    def __init__(self, modelFolder):
        self.modelDir = modelFolder

    @staticmethod
    def parseHeader(filename):
        symbols = []
        with open(filename, 'r') as f:
            for line in f:
                splitted = line.split(' ')
                if len(splitted) == 3:
                    symbols.append(splitted[1])
                else:
                    print('Unexpected line format in %s: %s' %
                          (filename, line))
        return symbols

    def readParameterNames(self):
        """Read list of parameter names"""
        return self.parseHeader(os.path.join(self.modelDir, 'p.h'))

    def readFixedParameterNames(self):
        """Read list of fixed parameter names"""
        return self.parseHeader(os.path.join(self.modelDir, 'k.h'))

    def readObservableNames(self):
        """Read list of observable names"""
        return self.parseHeader(os.path.join(self.modelDir, 'y.h'))

    def readObservableFormula(self):
        """Read list of observable formula"""
        import glob
        srcs = glob.glob(os.path.join(self.modelDir, '*_y.cpp'))
        assert len(srcs) == 1

        observables = []
        with open(srcs[0], 'r') as f:
            for line in f:
                match = re.sub(r'^ +y\[\d+\] = ([^;]+);\W$', r'\1', line, 1)
                if match == line:
                    continue
                observables.append(match)
        return observables

    def readObservables(self):
        # for backwards compatibility
        return self.readObservableFormula()


class HDF5DataGenerator:
    """
    Generate HDF5 file with fixed parameters and measurements for an AMICI-imported SBML model
    based on SBML model, AMICI syms file or model folder, measurement file and fixed parameter file.
    """

    def __init__(self, fileNameSBML, fileNameSyms, fileMeasurements, fileFixedParameters):
        """
        fileNameSBML: filename of model SBML file
        fileNameSyms: filename of AMICI model syms file or the model folder of a python-amici generated model
        fileMeasurements: filename of measurement file (TODO: specify format)
        fileFixedParameters: filename with AMICI fixed parameter vectors for all conditions refered to in the measurement file
        """

        # read symbols from compiled model to ensure correct ordering (might
        # change during sbml import)
        if fileNameSyms.endswith('.m'):
            # Matlab syms file, use AMICI model file interface for parsing
            self.amiciSyms = amiciHelper.AmiciSyms(fileNameSyms)
        else:
            # Python AMICI generated model folder
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

        print("Measurements shape", self.measurementDf.shape)

        # Get list of used conditions
        # Must also consider conditionRef values, since pre-equilibration conditions might come without any datapoints,
        # and thus, not show up in condition-column
        self.conditions = amiciHelper.unique(self.measurementDf.loc[:, 'condition'].append(
            self.measurementDf.loc[:, 'conditionRef']))
        # might have been introduced through empty conditionRef fields
        self.conditions = [c for c in self.conditions if not (
            isinstance(c, float) and math.isnan(c))]
        self.numConditions = len(self.conditions)

        # when using adjoint sensitivities, we cannot keep inf -> consider late
        # timepoint as steady-state
        print("Changing t = Inf to t = 1e8.")
        self.measurementDf.loc[self.measurementDf['time']
                               == np.inf, 'time'] = 1e8

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
        self.fixedParametersDf.set_index(
            self.fixedParametersDf.columns.values[0], inplace=True)

        print("Fixed parameters orginal: ", self.fixedParametersDf.shape)

        # drop conditions that do not have measurements
        dropCols = [
            label for label in self.fixedParametersDf if not label in self.conditions]
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

        self.generateReferenceConditionMap()

        print("Generating fixed parameters matrix...")
        self.generateFixedParameterMatrix()

        print("Generating measurement matrix...")
        self.generateMeasurementMatrix()

        print("Handling scaling parameters...")
        self.generateHierarchicalOptimizationData()

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
        self.writeStringArray(
            "/parameters/modelParameterNames", simulationParameterNames)
        print("Number of simulation parameters: %d" %
              len(simulationParameterNames))

        # simulationParameterMap is OrderedDict with simulation parameter names as keys, and dicts with
        # condition names as key and condition-specific instances as values.
        # no instance name means name is equal to model parameter name
        self.simulationParameterMap = collections.OrderedDict()
        for p in simulationParameterNames:
            self.simulationParameterMap[p] = {}

        # generate list of optimization parameters (>= simulation parameters
        # because of condition-specific parameters)
        self.optimizationParameterNames = self.replaceScalingParameterNames(
            simulationParameterNames[:])
        # keep map for reverse lookup
        self.optimizationParameterNamesToIndices = {
            name: idx for idx, name in enumerate(self.optimizationParameterNames)}

        print("Number of optimization parameters: %d" %
              len(self.optimizationParameterNames))
        self.writeStringArray("/parameters/parameterNames",
                              self.optimizationParameterNames)
        assert(len(simulationParameterNames) <=
               len(self.optimizationParameterNames))

        self.generateSimulationToOptimizationParameterMapping(
            simulationParameterNames, self.optimizationParameterNames)

        self.f.flush()

    def replaceScalingParameterNames(self, parameterNames):
        """
        Remove the "generic" names of the scaling parameter in the model and replace with the instances for the respective experiments
        """

        conditionSpecificScalingsUsed = self.getUsedScalingParameters()
        # print(scalingsUsed)
        for currentConditionSpecificName in conditionSpecificScalingsUsed:
            currentGenericName = self.getGenericParameterName(
                currentConditionSpecificName)
            # print(currentGenericName)
            if currentGenericName in parameterNames:
                parameterNames.remove(currentGenericName)

        parameterNames += conditionSpecificScalingsUsed
        return parameterNames

    def getUsedScalingParameters(self):
        """
        Get unique list of scaling parameters mentioned in measurement file
        with optimization parameter names (*not* as in sbml/amici model)
        """

        # List of condition-specific parameter names
        # NOTE: not using set() here which would scramble parameter order.
        # This allows use to keep starting points from /randomStarts
        # the same after regenerating the data file.
        conditionSpecificScalingParameterNames = []

        # Track if there have been condition-specific names provided for this
        # parameters. They need to be always or never condition-specific.
        observableHasConditionSpecificParameters = {}
        for index, row in self.measurementDf.iterrows():
            conditionName = row['condition']
            conditionSpecificNames = self.splitScalingParameterNames(
                row['scalingParameter'])
            if len(conditionSpecificNames):
                # ensure this observable has not been encountered or was marked
                # as conditonSpecific before
                if row['observable'] in observableHasConditionSpecificParameters and observableHasConditionSpecificParameters[row['observable']] == False:
                    raise ValueError(
                        "%s had condition-specific parameters before, but none are given for %s" % (row['observable'], row))
                else:
                    observableHasConditionSpecificParameters[row['observable']] = True
            else:
                if row['observable'] in observableHasConditionSpecificParameters and observableHasConditionSpecificParameters[row['observable']] == True:
                    raise ValueError(
                        "%s had no condition-specific parameters before, but has them in %s" % (row['observable'], row))
                else:
                    observableHasConditionSpecificParameters[row['observable']] = False

            for currentConditionSpecificName in conditionSpecificNames:
                conditionSpecificScalingParameterNames.append(
                    currentConditionSpecificName)
                currentGenericName = self.getGenericParameterName(
                    currentConditionSpecificName)
                self.simulationParameterMap[currentGenericName][conditionName] = currentConditionSpecificName

        conditionSpecificScalingParameterNames = amiciHelper.unique(
            conditionSpecificScalingParameterNames)
        return conditionSpecificScalingParameterNames

    def splitScalingParameterNames(self, names):
        if not isinstance(names, float):
            # N/A will be float
            # is single value or comma-separated list
            return names.split(",")
        return []

    def generateSimulationToOptimizationParameterMapping(self, simulationParameterNames, optimizationParameterNames):
        """
        Create dataset n_pararameters_simulation x n_conditions with indices of respective parameters in pararameters_optimization
        """

        numSimulationParameters = len(simulationParameterNames)

        # use in-memory matrix, don't write every entry to file directly (super
        # slow)
        parameterMap = np.zeros(
            shape=(numSimulationParameters, self.numConditions),)

        # create inverse mapping
        optimizationParameterNameToIndex = {
            name: idx for idx, name in enumerate(optimizationParameterNames)}

        for idxSimulationParameter, simulationParameterName in enumerate(simulationParameterNames):

            if not len(self.simulationParameterMap[simulationParameterName]):
                # not a condition-specific parameter, set same index for all
                # conditions
                parameterMap[idxSimulationParameter,
                             :] = optimizationParameterNameToIndex[simulationParameterName]
            else:
                # condition-specific parameter
                for conditionIdx, conditionName in enumerate(self.conditions):
                    # this must be present if there is data for the observable using this parameter (checked above),
                    # otherwise we can set any index, because amici will not use this parameter and set the gradient to 0.0. we will use 0 here.
                    # TODO: handle better?
                    if conditionName in self.simulationParameterMap[simulationParameterName]:
                        conditionSpecificParameterName = self.simulationParameterMap[
                            simulationParameterName][conditionName]
                        parameterMap[idxSimulationParameter,
                                     conditionIdx] = optimizationParameterNameToIndex[conditionSpecificParameterName]
                    else:
                        # This condition does not use the respective scaling parameter so there is no corresponding optimization parameter
                        # Cannot set to NaN in integer matrix. Will use 0.
                        # AMICI will not use this parameter anyways and its
                        # gradient will be 0.0
                        parameterMap[idxSimulationParameter, conditionIdx] = 0

        # write to file
        self.f.require_dataset('/parameters/optimizationSimulationMapping',
                               shape=(numSimulationParameters,
                                      self.numConditions),
                               chunks=(numSimulationParameters, 1),
                               dtype='<i4',
                               fillvalue=0,
                               compression=self.compression,
                               data=parameterMap)

    def getOptimizationParameterNameForConditionSpecificSimulationParameter(self, conditionIdx, simulationParameterName):
        """"For the given model parameter name and condition index, find the condition-specific optimization parameter"""
        # Measurements for given condition index
        scalingsForCurrentConditionByMeasurement = self.measurementDf.loc[
            self.measurementDf['condition'] == self.conditions[conditionIdx], 'scalingParameter']
        for scalingsForCurrentMeasurement in scalingsForCurrentConditionByMeasurement:
            for scaling in self.splitScalingParameterNames(scalingsForCurrentMeasurement):
                if self.getGenericParameterName(scaling) == simulationParameterName:
                    return scaling
        return None

    def getGenericParameterName(self, modelParameterName):
        """
        For the given condition-specific parameter name from measurement file or optimization parameter names,
        get the parameter name as specified in the sbml/amici model.

        Condition-specific parameter names are assumed to be named ${modelParameterName}_${suffixWithoutUnderscore},
        thus, remove suffix from given condition-specific parameter name.
        """

        return "_".join(modelParameterName.split("_")[:-1])

    def getGenericParameterNameAndCondition(self, scaling):
        """
        For the given condition-specific parameter name from measurement file or optimization parameter names,
        get the parameter name as specified in the sbml/amici model.

        Condition-specific parameter names are assumed to be named ${modelParameterName}_${suffixWithoutUnderscore},
        thus, split and return these two parts.
        """

        splitted = scaling.split("_")
        genericName = "_".join(splitted[:-1])
        suffix = splitted[-1]
        return genericName, suffix

    def getGenericParameterNames(self, conditionSpecificParameterNames):
        """
        For the given condition-specific parameter names from measurement file or optimization parameter names,
        return unique list of the parameter names as specified in the sbml/amici model.
        """
        # NOTE: not using set() here which would scramble parameter order.
        # This allows use to keep starting points from /randomStarts
        # the same after regenerating the data file.
        genericNames = []
        for s in conditionSpecificParameterNames:
            genericName = self.getGenericParameterName(s)
            genericNames.append(genericName)
        return amiciHelper.unique(genericNames)

    def generateFixedParameterMatrix(self):
        """
        Write fixed parameters dataset (nFixedParameters x nConditions).
        """

        k = self.amiciSyms.readFixedParameterNames()
        self.nk = len(k)
        print("Number of fixed parameters: %d" % len(k))

        # Create in-memory table, write all at once for speed
        fixedParameterMatrix = np.full(
            shape=(self.nk, self.numConditions), fill_value=np.nan)
        for i in range(len(k)):
            self.handleFixedParameter(i, k[i], fixedParameterMatrix)

        self.createFixedParameterDatasetAndWriteAttributes(
            k, fixedParameterMatrix)

        self.f.flush()

    def createFixedParameterDatasetAndWriteAttributes(self, fixedParameters, data):
        """
        Create fixed parameters data set and annotations
        """
        self.f.require_group("/fixedParameters")

        self.writeStringArray(
            "/fixedParameters/parameterNames", fixedParameters)
        self.writeStringArray(
            "/fixedParameters/conditionNames", self.conditions)

        # chunked for reading experiment-wise
        nk = len(fixedParameters)
        dset = self.f.create_dataset("/fixedParameters/k",
                                     (nk, self.numConditions),
                                     dtype='f8', chunks=(nk, 1), compression=self.compression, data=data)

        # set dimension scales
        dset.dims.create_scale(
            self.f['/fixedParameters/parameterNames'], 'parameterNames')
        dset.dims.create_scale(
            self.f['/fixedParameters/conditionNames'], 'conditionNames')
        dset.dims[0].attach_scale(self.f['/fixedParameters/parameterNames'])
        dset.dims[1].attach_scale(self.f['/fixedParameters/conditionNames'])

        return dset

    def generateReferenceConditionMap(self):
        """
        Write index map for reference conditions to be used for pre-equilibration

        Parse conditionRef column, and write to hdf5 dataset. If conditionRef is empty, set to -1, otherwise to condition index
        """

        # TODO: self.measurementDf.conditionRef: this should go to exp table to avoid potential redundancy/ambiguity in measurement file
        # -> really? could we want different preequilibration for a given condition?
        # TODO: we might want to separate reference (control, ...) condition
        # from preequilibration condition

        referenceMap = [NO_PREEQ_CONDITION_IDX] * len(self.conditions)

        for index, row in self.measurementDf.iterrows():
            conditionIdx = self.conditions.index(row['condition'])
            if isinstance(row['conditionRef'], float) and np.isnan(row['conditionRef']):
                conditionIdxRef = NO_PREEQ_CONDITION_IDX
            else:
                conditionIdxRef = self.conditions.index(row['conditionRef'])

            if referenceMap[conditionIdx] != conditionIdxRef and referenceMap[conditionIdx] != NO_PREEQ_CONDITION_IDX:
                print("Error: Ambiguous assignment of reference conditions for %s: %s vs. %s"
                      % (self.conditions[conditionIdx], self.conditions[referenceMap[conditionIdx]], self.conditions[conditionIdxRef]))

            referenceMap[conditionIdx] = conditionIdxRef

        self.f.create_dataset("/fixedParameters/referenceConditionsX",
                              shape=(self.numConditions,), dtype="<i4", data=referenceMap)

    def writeStringArray(self, path, strings):
        """
        Write string array to hdf5
        """
        dt = h5py.special_dtype(vlen=str)
        dset = self.f.create_dataset(path, (len(strings),), dtype=dt)
        dset[:] = [s.encode('utf8') for s in strings]
        self.f.flush()

    def writeFloatArray(self, path, values, dtype='f8'):
        """
        Write float array to hdf5
        """
        dset = self.f.create_dataset(path, (len(values),), dtype=dtype)
        dset[:] = values
        self.f.flush()

    def writeIntArray(self, path, values, dtype='<i4'):
        """
        Write integer array to hdf5
        """
        dset = self.f.create_dataset(path, (len(values),), dtype=dtype)
        dset[:] = values
        self.f.flush()

    def handleFixedParameter(self, parameterIndex, parameterName, dset):
        """
        Extract parameter values from data table or model and set to dset

        Lookup order is condition table > model parameter > model species

        condition table is expected to have "globalized names", i.e. reactionId_localParameterId

        NOTE: parameters which are not provided in fixedParametersDf or model are set to 0.0

        Arguments:
        parameterIndex: Index of the current fixed parameter
        parameterName: Name of the current parameter
        dset: 2D array-like (nFixedParameters x nConditions) where the parameter value is to be set
        """

        if parameterName in self.fixedParametersDf.index:
            # Parameter in condition table
            dset[parameterIndex, :] = self.fixedParametersDf.loc[parameterName,
                                                                 self.conditions].values
        else:
            sbml_parameter = self.sbmlModel.getParameter(parameterName)
            if sbml_parameter:
                # Parameter value from model
                dset[parameterIndex, :] = sbml_parameter.getValue()
            else:
                sbml_species = self.sbmlModel.getSpecies(parameterName)
                if sbml_species:
                    # A constant species might have been turned in to a model parameter
                    # TODO: we dont do any conversion here, although we would want to have concentration
                    # currently there is only 1.0
                    dset[parameterIndex, :] = sbml_species.getInitialConcentration() if sbml_species.isSetInitialConcentration() else sbml_species.getInitialAmount()
                else:
                    # We need to check for "globalized" parameter names too (reactionId_localParameterId)
                    # model has localParameterId, data file has globalized name
                    global_name = self.getGlobalNameForLocalParameter(parameterName)
                    if global_name:
                        sbml_parameter = self.sbmlModel.getParameter(global_name)
                        if sbml_parameter:
                            # Use model parameter value
                            dset[parameterIndex, :] = sbml_parameter.getValue()
                        else:
                            print("Warning: Fixed parameter not found in ExpTable, setting to 0.0: ", parameterName)
                            dset[parameterIndex, :] = 0.0
                    else:
                        print("Warning: Fixed parameter not found in ExpTable, setting to 0.0: ", parameterName)
                        dset[parameterIndex, :] = 0.0

    def getGlobalNameForLocalParameter(self, needle_parameter_id):
        sbml_model = self.sbmlModel
        for reaction in sbml_model.getListOfReactions():
            kl = reaction.getKineticLaw()
            for p in kl.getListOfParameters():
                parameter_id = p.getId()
                if parameter_id.endswith(needle_parameter_id):
                        return f'{reaction.getId()}_{parameter_id}'
        return None


    def generateMeasurementMatrix(self):
        """
        Generate matrix with training data for all observables defined in the _syms file.

        NOTE: data is expected to have numTimepoints as leading dimensions. parameter estimation might fail silently if not.
        """

        self.f.create_group("/measurements")

        self.y = self.amiciSyms.readObservables()
        self.ny = len(self.y)
        self.observables = self.getObservableNamesFromSbml()
        self.checkObservablesSbmlMatchAmici(self.observables, self.y)
        self.writeStringArray(
            "/measurements/observableNames", self.observables)

        print("Number of observables: %d" % self.ny)

        dsetY = self.f.create_dataset("/measurements/y",
                                      shape=(self.numConditions,
                                             self.numTimepoints, self.ny),
                                      chunks=(1, self.numTimepoints, self.ny),
                                      fillvalue=np.nan, dtype='f8', compression=self.compression)

        dsetSigmaY = self.f.create_dataset("/measurements/ysigma",
                                           shape=(self.numConditions,
                                                  self.numTimepoints, self.ny),
                                           chunks=(
                                               1, self.numTimepoints, self.ny),
                                           fillvalue=np.nan, dtype='f8', compression=self.compression)

        self.writeMeasurements(dsetY, dsetSigmaY)
        self.f.flush()

        # TODO: this will be addressed when enabling preequilibration, for now there is the referenceConditionMap
        # one is measurement reference -> need for generating simulation data
        # dset = self.f.create_dataset("/fixedParameters/reference",
        #                             (nk, self.numConditions),
        #                             dtype='f8', chunks=(nk, 1))
        # reference = map[condition][observable]
        # self.writeStringArray("/fixedParameters/celllineNames", self.cells)
        # self.writeIntArray("/fixedParameters/celllineIdx", self.getConditionCelllineIndexMap())
        # g.attrs['numCelllines'] = len(set(self.cells))

    def getObservableNamesFromSbml(self):
        """
        Get array with names of observables from SBML model.

        NOTE: Observables are not accounted for in the SBML standard.
        We implement them as parameters and assignment rules.
        The respective parameters start with "observable_",
        but do not end in "_sigma", which is reserved for encoding the error model.

        NOTE: observable names in the model are expected to match the ones in the data table
        """
        observables = []
        for p in self.sbmlModel.getListOfParameters():
            p = p.getId()
            if p.startswith('observable_') and not p.endswith('_sigma'):
                #observableName = p[len('observable_'):]
                observableName = p
                observables.append(observableName)
        observablesUni = amiciHelper.unique(observables)

        if len(observables) != len(observablesUni):
            print("!! %d redundant observable definitions in model"
                  % (len(observables) - len(observablesUni)))
        return observablesUni

    def writeMeasurements(self, dsetY, dsetSigmaY):
        """
        Write measurements to hdf5 dataset
        """
        # TODO: at the moment we cannot deal with replicate measurements (same condition, same observable, same timepoint) since not supported by AMICI
        # should alert user in case that happens, currently will be silently overwritten by the latter
        # dirty workaround: take mean of change time to t+eps

        # write data from each row in measurementDf
        for index, row in self.measurementDf.iterrows():
            conditionIdx = self.conditions.index(row['condition'])
            observableIdx = self.observables.index(row['observable'])
            timeIdx = self.timepoints.index(row['time'])
            dsetY[conditionIdx, timeIdx, observableIdx] = float(
                row['measurement'])

            scalings = self.splitScalingParameterNames(row['scalingParameter'])
            hasSigmaParameter = len(
                [s for s in scalings if self.getGenericParameterName(s).endswith('_sigma')])
            sigma = float(row['sigma'])
            if hasSigmaParameter and not np.isnan(sigma):
                print(colored('Warning: Sigma parameter name and fixed sigma provided at the same time for "%s". Setting to NaN so that sigma parameter is used' % (row), 'yellow'))
                sigma = np.nan
            dsetSigmaY[conditionIdx, timeIdx, observableIdx] = sigma

    def checkObservablesSbmlMatchAmici(self, observables, y):
        """
        Check if order of observables in SBML model match that of AMICI model (which is assumed here)

        Problem: order may have changed due to SBML which does not necessarily conserve ordering
        Problem: formula matching is not straightforward, since they might have changed due to symbolic processing

        TODO: require finished amici python module for data generation, get order from there instead of sbml model
        This will get rid of the checks
        """

        # The least we can do is compare the number of observables
        if len(y) != len(observables):
            raise AssertionError('y: %d, obs: %d' % (len(y), len(observables)))

        for i in range(len(observables)):
            observableName = observables[i]
            amiciObsFormula = y[i]

            # rough check if the sbml and model observable match
            # this check expected observables to be names
            # observable_$observedSpecies
            obsSpecies = observableName[(observableName.rfind('_') + 1):]
            if amiciObsFormula.find(obsSpecies) < 0:
                print(colored('Warning: Not sure if SBML observable "%s" matches AMICI formula "%s"' % (
                    observableName, amiciObsFormula), 'yellow'))

    def generateHierarchicalOptimizationData(self):
        """
        Deal with offsets, proportionality factors and sigmas for hierarchical optimization

        Generate the respective index lists and mappings
        """

        # TODO: should check at the end of the function if scalingIndices lists
        # are non-overlapping
        self.ensureNonOverlappingParameterForHierarchicalOptimization()

        self.handleProportionalityFactors()
        self.handleOffsetParameter()  # must call after handleProportionalityFactors
        self.handleSigmas()

    def ensureNonOverlappingParameterForHierarchicalOptimization(self):
        """
        Current parPE implementation of hierarchical optimization assumes that only one single
        offset or proportionality parameter is to be calculated per condition x observable.
        Will fail silently if there are more. Thus we need to check here.
        Also check that there is at most one sigma per condition x observable
        """

        parameterMap = {}  # observable => condition => parameter-for-hierarchical
        for index, row in self.measurementDf.iterrows():
            scalings = self.splitScalingParameterNames(row['scalingParameter'])
            if not len(scalings):
                continue

            nonSigmas = [x for x in scalings if not self.getGenericParameterName(
                x).endswith("_sigma")]
            sigmas = [x for x in scalings if self.getGenericParameterName(
                x).endswith("_sigma")]

            if len(sigmas) > 1:
                raise RuntimeError(
                    'ERROR: Multiple sigma parameters provided for\n%s\n' % row)

            if len(nonSigmas) > 1:
                print("Warning: multiple scaling and/or offset parameters found for\n\t%s\n\twhich one to choose for hierarchical optimization? Assuming first" % row)
            elif not len(nonSigmas):
                continue

            if row['observable'] in parameterMap:
                if row['condition'] in row['observable'] and parameterMap[row['observable']][row['condition']] != nonSigmas[0]:
                    raise RuntimeError("ERROR: multiple scaling and/or offset parameters found for\n\t%s\n\twhich one to choose for hierarchical optimization? previously found %s" % (
                        row, parameterMap[row['observable']]))
            else:
                parameterMap[row['observable']] = {}
            parameterMap[row['observable']][row['condition']] = nonSigmas[0]

    def handleOffsetParameter(self):
        """
        Write list of offset parameters selected for hierarchical optimization
        """

        '''
        # list of observables for which there are already analytically computed parameters
        observablesBlacklist = set()
        if '/scalingParametersMapToObservables' in self.f:
            observablesBlacklist = set(self.f["/scalingParametersMapToObservables"][:, 2])
                                      
        # TODO ensure that this observable does not have a proportionality factor
        #offsetsForHierarchical = [x for x in offsetsForHierarchical if x.startswith("offset_") and getObservableForScalingParameter not in observablesBlacklist ]
        '''
        offsetsForHierarchical = [
            x for x in self.getUsedScalingParameters() if x.startswith("offset_")]
        # don't create dataset if it would be empty
        if not len(offsetsForHierarchical):
            return

        # find indices for names
        offsetsForHierarchicalIndices = [
            self.optimizationParameterNamesToIndices[x] for x in offsetsForHierarchical]

        [self.ensureOffsetIsOffsetWithPositiveSign(
            x) for x in offsetsForHierarchical]
        print("Number of offset parameters for hierarchical optimization: %d" %
              len(offsetsForHierarchicalIndices))

        self.f.require_dataset("/offsetParameterIndices",
                               shape=(len(offsetsForHierarchicalIndices),),
                               dtype='<i4',
                               data=offsetsForHierarchicalIndices)

        # find usages for the selected parameters
        use = self.getAnalyticalParameterTable(offsetsForHierarchical)

        self.f.require_dataset("/offsetParametersMapToObservables",
                               shape=(len(use), 3),
                               dtype='<i4', data=use)

    def getAnalyticalParameterTable(self, parametersForHierarchical):
        """Generate (scalingIdx, conditionIdx, observableIdx) table for all occurrences of the given parameter names

        Returns:
        list of (scalingIdx, conditionIdx, observableIdx) tuples
        """
        use = []
        for index, row in self.measurementDf.iterrows():
            currentScalings = self.splitScalingParameterNames(
                row['scalingParameter'])
            for s in currentScalings:
                #print(s, parametersForHierarchical)
                try:
                    scalingIdx = parametersForHierarchical.index(s)
                    conditionIdx = self.conditions.index(row['condition'])
                    observableIdx = self.observables.index(row['observable'])
                    tup = (scalingIdx, conditionIdx, observableIdx)
                    # Don't add a new line for each timepoint
                    # We don't allow separate parameters for individual time-points
                    # (Can be implemented via different observables)
                    if not tup in use:
                        use.append(tup)
                except ValueError:
                    pass  # current parameter not in list
        assert len(use)
        return use

    def ensureOffsetIsOffsetWithPositiveSign(self, scaling):
        """
        Current parPE implementation of hierarchical optimization assumes that offset parameters have positive sign.
        Will fail silently if this is not true, therefore check here that it is in indeed an offset parameter and
        that it has positive sign.
        """

        name, formulae = self.getNameAndFormulasForConditionSpecificParameter(
            scaling)
        param = sp.sympify(name)
        itsOkay = True
        for formula in formulae:
            obs = sp.sympify(formula)
            itsOkay = itsOkay and (param not in (obs - param).free_symbols)

            if itsOkay:
                print(colored("Parameter %s selected as offset for hierarchical optimization (%s)." % (
                    scaling, formula), "green"))
            else:
                print(colored("ERROR: Parameter %s incorrectly selected as offset factor for hierarchical optimization (%s). Has positive sign?" % (
                    scaling, formula), "green"))

        if not itsOkay:
            raise AssertionError

    def getNameAndFormulasForConditionSpecificParameter(self, parameterName):
        """
        Get formula for the (first) observable associated with the provided condition-specific parameter
        (name provided as in measurement list, not as in model)
        """
        # get model parameter name
        name = self.getGenericParameterName(parameterName)
        formulaMatches = []
        for formula in self.y:
            if re.search(r'(^|\W)' + re.escape(name) + r'($|\W)', formula):
                formulaMatches.append(formula)
        return name, formulaMatches

    def handleProportionalityFactors(self):
        """
        Write datasets specifying which proportionality factors to consider for hierarchical optimization
        """

        # TODO: Figure out how to deal with observables with multiple scaling
        # parameters; which one is preferred? proportionality over offset? how
        # to let the user select which ones to add?

        scalingsForHierarchical = [
            x for x in self.getUsedScalingParameters() if x.startswith("scaling_")]
        if not len(scalingsForHierarchical):
            return
        scalingsForHierarchicalIndices = [
            self.optimizationParameterNamesToIndices[x] for x in scalingsForHierarchical]
        [self.ensureProportionalityFactorIsProportionalityFactor(
            x) for x in scalingsForHierarchical]

        self.f.require_dataset("/scalingParameterIndices",
                               shape=(len(scalingsForHierarchicalIndices),),
                               dtype='<i4',
                               data=scalingsForHierarchicalIndices)
        print("Number of proportionality factors for hierarchical optimization: %d" % len(
            scalingsForHierarchicalIndices))

        # find usages for the selected parameters
        use = self.getAnalyticalParameterTable(scalingsForHierarchical)

        self.f.require_dataset("/scalingParametersMapToObservables",
                               shape=(len(use), 3),
                               dtype='<i4', data=use)

    def ensureProportionalityFactorIsProportionalityFactor(self, scaling):
        """
        Ensure that this is a proportionality factor (a as in y = a*x)
        Note: a*x + b currently not supported
        """

        name, formulae = self.getNameAndFormulasForConditionSpecificParameter(
            scaling)
        param = sp.sympify(name)
        itsOkay = True
        for formula in formulae:
            obs = sp.sympify(formula)
            itsOkay = itsOkay and (param not in (obs / param).free_symbols)

            if itsOkay:
                print(colored("Parameter %s selected as proportionality factor for hierarchical optimization (%s)." % (
                    scaling, formula), "green"))
            else:
                print(colored("ERROR: Parameter %s incorrectly selected as proportionality factor for hierarchical optimization (%s)." % (
                    scaling, formula), "green"))

        if not itsOkay:
            raise AssertionError

    def handleSigmas(self):
        """
        Write data for dealing with sigma parameters in hierarchical optimization
        """
        sigmasForHierarchical = [x for x in self.getUsedScalingParameters(
        ) if self.getGenericParameterName(x).endswith("_sigma") or self.getGenericParameterName(x).startswith("sigma_")]
        if not len(sigmasForHierarchical):
            return

        sigmasForHierarchicalIndices = [
            self.optimizationParameterNamesToIndices[x] for x in sigmasForHierarchical]

        self.f.require_dataset("/sigmaParameterIndices",
                               shape=(len(sigmasForHierarchicalIndices),),
                               dtype='<i4',
                               data=sigmasForHierarchicalIndices)
        print("Number of sigmas for hierarchical optimization: %d" %
              len(sigmasForHierarchicalIndices))

        # find usages for the selected parameters
        use = self.getAnalyticalParameterTable(sigmasForHierarchical)

        self.f.require_dataset("/sigmaParametersMapToObservables",
                               shape=(len(use), 3),
                               dtype='<i4', data=use)

    def copyAmiciOptions(self):
        """
        Write simulation options
        """
        g = self.f.require_group("/amiciOptions")
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

        numParameters = self.f['/parameters/modelParameterNames'].shape[0]

        # parameter indices w.r.t. which to compute sensitivities
        self.f.require_dataset(
            '/amiciOptions/sens_ind', shape=(numParameters,), dtype="<i4", data=range(numParameters))

        self.f.require_dataset(
            '/amiciOptions/ts', shape=(len(self.timepoints),), dtype="f8", data=self.timepoints)

        # set parameter scaling: all log10, except for offsets which can be negative
        # ... for AMICI model parameters:
        offsetIndices = [i for i, p in enumerate(
            self.f['/parameters/modelParameterNames']) if p.startswith('offset_')]
        pscale = np.full(numParameters, 2)
        pscale[offsetIndices] = 0
        self.f.require_dataset('/amiciOptions/pscale',
                               shape=(numParameters,), dtype="<i4", data=pscale)

        # ... for all parameters for hierarchical optimization
        offsetIndices = [i for i, p in enumerate(
            self.f['/parameters/parameterNames']) if p.startswith('offset_')]
        # self.getAnalyticallyComputedSimulationParameterIndices()
        linParametersAmiciIndices = offsetIndices
        numOptimizationParameters = self.f['/parameters/parameterNames'].shape[0]
        self.f.require_dataset('/parameters/pscale',
                               shape=(numOptimizationParameters,), dtype="<i4", data=[2 * (ip not in linParametersAmiciIndices) for ip in range(numOptimizationParameters)])

    def getAnalyticallyComputedSimulationParameterIndices(self):
        """
        Get model parameter index (not optimization parameter index) of all analytically computed parameters
        """
        parameterNamesModel = []
        if '/offsetParameterIndices' in self.f:
            parameterNamesOptimization = self.f['/parameters/parameterNames'][self.f['/offsetParameterIndices']]
            parameterNamesModel.extend(
                set([self.getGenericParameterName(o) for o in parameterNamesOptimization]))

        if '/scalingParameterIndices' in self.f:
            parameterNamesOptimization = self.f['/parameters/parameterNames'][self.f['/scalingParameterIndices']]
            parameterNamesModel.extend(
                set([self.getGenericParameterName(o) for o in parameterNamesOptimization]))

        if '/sigmaParameterIndices' in self.f:
            parameterNamesOptimization = self.f['/parameters/parameterNames'][self.f['/sigmaParameterIndices']]
            parameterNamesModel.extend(
                set([self.getGenericParameterName(o) for o in parameterNamesOptimization]))

        return [self.f['/parameters/modelParameterNames'][:].tolist().index(p) for p in set(parameterNamesModel)]

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
        g.attrs['optimizer'] = 0  # IpOpt
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
        # set ridiculously high, so only the acceptable_* options below matter
        g.attrs["acceptable_tol"] = 1e20
        g.attrs["acceptable_obj_change_tol"] = 1e-12
        g.attrs["watchdog_shortened_iter_trigger"] = 0

        # set fmincon options
        g = self.f.require_group('optimizationOptions/fmincon')
        g.attrs['MaxIter'] = 100
        g.attrs["TolX"] = 1e-8
        g.attrs["TolFun"] = 0
        g.attrs["MaxFunEvals"] = 1e7
        g.attrs["algorithm"] = np.string_("interior-point")
        g.attrs["GradObj"] = np.string_("on")
        g.attrs["display"] = np.string_("iter")

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

        Offset parameters are allowed to be negative
        """
        numParams = self.f['/parameters/parameterNames'].shape[0]
        lower = self.f.require_dataset(
            '/parameters/lowerBound', [numParams], 'f8')
        upper = self.f.require_dataset(
            '/parameters/upperBound', [numParams], 'f8')
        lower[:] = [-5] * numParams
        upper[:] = [3] * numParams

        # offset parameters are optimized in linear space
        offsetIndices = [i for i, p in enumerate(
            self.f['/parameters/parameterNames']) if p.startswith('offset_')]
        if len(offsetIndices):
            lower[offsetIndices] = -1e10
            upper[offsetIndices] = +1e10

    def writeStartingPoints(self):
        """
        Write a list of random starting points uniformly sampled from the parameter bounds.
        Parameter bounds need to be written beforehand.
        """
        numParams = self.f['/parameters/parameterNames'].shape[0]
        numStartingPoints = 100
        np.random.seed(0)
        startingPoints = self.f.require_dataset(
            '/optimizationOptions/randomStarts', [numParams, numStartingPoints], 'f8')
        lower = self.f['/parameters/lowerBound'][:]
        upper = self.f['/parameters/upperBound'][:]
        startingPoints[:] = np.transpose(np.random.rand(numStartingPoints, numParams) * (upper - lower) + lower)


def main():
    if len(sys.argv) < 6:
        print("Usage: %s hdf5File sbmlModelFile symsModelFile/pyModelFolder cost_func exp_table" % __file__)
        sys.exit()

    (hdf5File, sbmlModelFile, symsModelFile,
     costFunctionFile, expTableFile) = sys.argv[1:6]
    h5gen = HDF5DataGenerator(
        sbmlModelFile, symsModelFile, costFunctionFile, expTableFile)
    h5gen.generateFile(hdf5File)


if __name__ == "__main__":
    main()
