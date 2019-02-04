#!/usr/bin/env python3
"""
Generate HDF5 file for parPE with fixed parameters and measurements for an
AMICI-imported SBML model based on tables with fixed parameters and training
data in the PEtab format https://github.com/ICB-DCM/PEtab

2018 Daniel Weindl <daniel.weindl@helmholtz-muenchen.de>
     Leonard Schmiester <leonard.schmiester@helmholtz-muenchen.de>

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
import argparse
import petab
import importlib

SCALING_LOG10 = 2
SCALING_LIN = 0
# index for no reference/preequilibration condition
NO_PREEQ_CONDITION_IDX = -1


def writeStringArray(f, path, strings):
    """
    Write string array to hdf5
    """
    dt = h5py.special_dtype(vlen=str)
    dset = f.create_dataset(path, (len(strings),), dtype=dt)
    dset[:] = [s.encode('utf8') for s in strings]
    f.flush()


def writeFloatArray(f, path, values, dtype='f8'):
    """
    Write float array to hdf5
    """
    dset = f.create_dataset(path, (len(values),), dtype=dtype)
    dset[:] = values
    f.flush()


def writeIntArray(f, path, values, dtype='<i4'):
    """
    Write integer array to hdf5
    """
    dset = f.create_dataset(path, (len(values),), dtype=dtype)
    dset[:] = values
    f.flush()


class HDF5DataGenerator:
    """
    Generate HDF5 file with fixed parameters and measurements for an
    AMICI-imported SBML model based on SBML model file, AMICI python model,
    measurement table, fixed parameter table and parameter table
    in PEtab format.

    Attributes:
        compression: h5py compression to be used
        measurement_df: pandas.DataFrame
            PEtab measurement table
        condition_df: pandas.DataFrame
            PEtab condition table
        parameter_df: pandas.DataFrame
            PEtab parameter table
        sbml_model:
        condition_ids: numpy.array condition IDs
        num_conditions: number of conditions
        timepoints: time points for which there is data
        num_timepoints: number of time points for which there is data
        f: h5py.File
            hdf5 file which is being created
    """

    def __init__(self, fileNameSBML, fileMeasurements, fileConditions,
                 fileParameters, model_output_dir, model_name):
        """
        fileNameSBML: filename of model SBML file (PEtab-style)
        fileMeasurements: filename of measurement file
        fileFixedParameters: filename with AMICI fixed parameter vectors for
            all conditions referred to in the measurement file
        fileParameters: PEtab parameter table filename
        """
        self.petab_problem = petab.Problem(
            sbml_file=fileNameSBML,
            condition_file=fileConditions,
            measurement_file=fileMeasurements,
            parameter_file=fileParameters
        )
        petab.lint_problem(self.petab_problem)

        # load input data
        self.parseMeasurementFile(fileMeasurements)
        self.parseFixedParametersFile(fileConditions)
        self.parameter_df = petab.get_parameter_df(fileParameters)
        self.loadSBMLModel(fileNameSBML)

        # load model
        sys.path.insert(0, model_output_dir)
        model_module = importlib.import_module(model_name)
        self.model = model_module.getModel()

        # scriptDir = os.path.dirname(os.path.realpath(__file__))

        # hdf5 dataset compression
        self.compression = "gzip"

    def parseMeasurementFile(self, filename):
        """
        Read cost_fun file and determine number of conditions and timepoints
        """
        self.measurement_df = petab.get_measurement_df(filename)

        print("Measurements shape:", self.measurement_df.shape)

        # Get list of used conditions
        self.condition_ids = amiciHelper.unique(
            self.measurement_df.simulationConditionId.values)
        self.num_conditions = len(self.condition_ids)

        # when using adjoint sensitivities, we cannot keep inf
        # -> consider late timepoint as steady-state
        print("Changing t = Inf to t = 1e8.")
        self.measurement_df.loc[
            self.measurement_df.time == np.inf, 'time'] = 1e8

        # list of timepoints
        # TODO: this will have to allowed to contain duplicates, if we want to
        # deal with replicate measurements
        self.timepoints = sorted(self.measurement_df.time.unique())
        self.num_timepoints = len(self.timepoints)

        print("Num conditions: ", self.num_conditions)
        print("Num timepoints: ", self.num_timepoints, self.timepoints)

    def parseFixedParametersFile(self, filename):
        """
        Load file and select relevant conditions
        """

        self.condition_df = pd.read_csv(filename, delimiter="\t")
        print("Fixed parameters orginal: ", self.condition_df.shape)

        # drop conditions that do not have measurements
        dropCols = [
            label for label in self.condition_df if
            not label in self.condition_ids]
        self.condition_df.drop(dropCols, 1, inplace=True)

        print("Fixed parameters usable: ", self.condition_df.shape)

    def loadSBMLModel(self, filename):
        """
        Load SBML model
        """
        # Note: must keep reader and document, otherwise model becomes invalid
        self.sbml_reader = SBMLReader()
        self.sbml_document = self.sbml_reader.readSBMLFromFile(filename)
        self.sbml_model = self.sbml_document.getModel()

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
        model_parameter_ids = np.array(self.model.getParameterIds())
        writeStringArray(self.f, "/parameters/modelParameterNames",
                         model_parameter_ids)
        print("Number of model parameters:", len(model_parameter_ids))

        print("Number of optimization parameters: %d" %
              len(self.parameter_df))
        writeStringArray(self.f, "/parameters/parameterNames",
                         self.parameter_df.index.name)
        # assert (len(model_parameter_ids) <=
        #        len(self.parameter_df))

        self.generateSimulationToOptimizationParameterMapping()

        self.f.flush()

    def generateSimulationToOptimizationParameterMapping(self):
        """
        Create dataset n_parameters_simulation x n_conditions with indices of
        respective parameters in pararameters_optimization
        """

        # get list of list of parameters for condition
        # outer has length of number of conditions
        # inner has length of model parameters
        self.parameter_mapping = self.petab_problem \
            .get_optimization_to_simulation_parameter_mapping()

        num_model_parameters = self.model.np()

        # use in-memory matrix, don't write every entry to file directly (super
        # slow)
        mapping_matrix = np.zeros(
            shape=(num_model_parameters, self.num_conditions), )

        # create inverse mapping for faster lookup
        optimization_parameter_name_to_index = {
            name: idx for idx, name in enumerate(self.parameter_df.index)}
        print(optimization_parameter_name_to_index)
        self.optimization_parameter_name_to_index = optimization_parameter_name_to_index

        model_parameter_ids = self.model.getParameterIds()

        for condition_idx, condition_map in enumerate(self.parameter_mapping):
            for model_parameter_idx, mapped_parameter in enumerate(
                    condition_map):
                # print(condition_idx, model_parameter_idx, mapped_parameter)
                if isinstance(mapped_parameter, str):
                    mapping_matrix[model_parameter_idx, condition_idx] = \
                        optimization_parameter_name_to_index[mapped_parameter]
                else:
                    # Have numeric override
                    # This condition does not use any parameter override here
                    # Cannot set to NaN in integer matrix. Will use 0.
                    # AMICI will not use this parameter anyways and its
                    # gradient will be 0.0
                    # numeric overrides for sigmas will be handled by model
                    # others are not yet supported
                    model_parameter_id = model_parameter_ids[
                        model_parameter_idx]
                    mapping_matrix[model_parameter_idx, condition_idx] = 0
                    print(f"Numeric override for condition {condition_idx} "
                          f"parameter {model_parameter_id}: "
                          f"{mapped_parameter}")
                    if not model_parameter_id.startswith('noiseParameter1'):
                        raise RuntimeError(
                            'Cannot handle numeric non-noise parameter overrides.')

        # write to file
        self.f.require_dataset(
            name='/parameters/optimizationSimulationMapping',
            shape=(num_model_parameters, self.num_conditions),
            chunks=(num_model_parameters, 1),
            dtype='<i4',
            fillvalue=0,
            compression=self.compression,
            data=mapping_matrix)

    def generateFixedParameterMatrix(self):
        """
        Write fixed parameters dataset (nFixedParameters x nConditions).
        """

        fixed_parameter_ids = self.model.getFixedParameterIds()
        self.nk = len(fixed_parameter_ids)
        print("Number of fixed parameters: %d" % len(fixed_parameter_ids))

        # Create in-memory table, write all at once for speed
        fixedParameterMatrix = np.full(
            shape=(self.nk, self.num_conditions), fill_value=np.nan)
        for i in range(len(fixed_parameter_ids)):
            self.handleFixedParameter(i, fixed_parameter_ids[i],
                                      fixedParameterMatrix)

        self.createFixedParameterDatasetAndWriteAttributes(
            fixed_parameter_ids, fixedParameterMatrix)

        self.f.flush()

    def createFixedParameterDatasetAndWriteAttributes(self,
                                                      fixed_parameter_ids,
                                                      data):
        """
        Create fixed parameters data set and annotations
        """
        self.f.require_group("/fixedParameters")

        writeStringArray(self.f, "/fixedParameters/parameterNames",
                         fixed_parameter_ids)
        writeStringArray(self.f,
                         "/fixedParameters/conditionNames", self.condition_ids)

        # chunked for reading experiment-wise
        nk = len(fixed_parameter_ids)
        dset = self.f.create_dataset("/fixedParameters/k",
                                     (nk, self.num_conditions),
                                     dtype='f8', chunks=(nk, 1),
                                     compression=self.compression, data=data)

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
        Write index map for conditions to be used for pre-equilibration.

        Parse preequilibrationConditionId column, and write to hdf5 dataset.
        If conditionRef is empty, set to -1, otherwise to condition index.

        """
        referenceMap = [NO_PREEQ_CONDITION_IDX] * self.num_conditions

        if np.issubdtype(self.measurement_df.preequilibrationConditionId.dtype,
                         np.number):
            # all NaNs, so there is no preeq required:
            pass
        else:
            # TODO: This assumes a one to one mapping of simulationCondition and
            # (possibly empty) preequilibrationCondition
            # we might have t allow for different combinations there, which will
            # extend the simulation condition vector over the number of unique
            # conditions listed in the measurement file
            assert len(self.measurement_df.simulationConditionId.unique()) \
                   == len(self.measurement_df.groupby(['simulationConditionId',
                                                       'preequilibrationConditionId']).size())

            for index, row in self.measurement_df.groupby(
                    ['simulationConditionId',
                     'preequilibrationConditionId']).size().reset_index().iterrows():
                conditionIdx = self.condition_ids.index(
                    row['simulationConditionId'])
                if isinstance(row['preequilibrationConditionId'],
                              float) and np.isnan(
                        row['preequilibrationConditionId']):
                    conditionIdxRef = NO_PREEQ_CONDITION_IDX
                else:
                    conditionIdxRef = self.condition_ids.index(
                        row['preequilibrationConditionId'])

                referenceMap[conditionIdx] = conditionIdxRef

        self.f.create_dataset("/fixedParameters/referenceConditionsX",
                              shape=(self.num_conditions,), dtype="<i4",
                              data=referenceMap)

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

        if parameterName in self.condition_df.index:
            # Parameter in condition table
            dset[parameterIndex, :] = self.condition_df.loc[parameterName,
                                                            self.condition_ids].values
        else:
            sbml_parameter = self.sbml_model.getParameter(parameterName)
            if sbml_parameter:
                # Parameter value from model
                dset[parameterIndex, :] = sbml_parameter.getValue()
            else:
                sbml_species = self.sbml_model.getSpecies(parameterName)
                if sbml_species:
                    # A constant species might have been turned in to a model parameter
                    # TODO: we dont do any conversion here, although we would want to have concentration
                    # currently there is only 1.0
                    dset[parameterIndex,
                    :] = sbml_species.getInitialConcentration() if sbml_species.isSetInitialConcentration() else sbml_species.getInitialAmount()
                else:
                    # We need to check for "globalized" parameter names too (reactionId_localParameterId)
                    # model has localParameterId, data file has globalized name
                    global_name = self.getGlobalNameForLocalParameter(
                        parameterName)
                    if global_name:
                        sbml_parameter = self.sbml_model.getParameter(
                            global_name)
                        if sbml_parameter:
                            # Use model parameter value
                            dset[parameterIndex, :] = sbml_parameter.getValue()
                        else:
                            print(
                                "Warning: Fixed parameter not found in ExpTable, setting to 0.0: ",
                                parameterName)
                            dset[parameterIndex, :] = 0.0
                    else:
                        print(
                            "Warning: Fixed parameter not found in ExpTable, setting to 0.0: ",
                            parameterName)
                        dset[parameterIndex, :] = 0.0

    def getGlobalNameForLocalParameter(self, needle_parameter_id):
        sbml_model = self.sbml_model
        for reaction in sbml_model.getListOfReactions():
            kl = reaction.getKineticLaw()
            for p in kl.getListOfParameters():
                parameter_id = p.getId()
                if parameter_id.endswith(needle_parameter_id):
                    return f'{reaction.getId()}_{parameter_id}'
        return None

    def generateMeasurementMatrix(self):
        """
        Generate matrix with training data for all observables and timepoints

        dim: num_conditions x num_timepoints x num_observables
        """

        self.f.create_group("/measurements")

        self.observable_ids = self.model.getObservableIds()
        # trim observable_ TODO: should be done in amici import
        self.observable_ids = [o[len('observable_'):] for o in
                               self.observable_ids]
        self.ny = self.model.ny
        writeStringArray(self.f,
                         "/measurements/observableNames", self.observable_ids)

        print("Number of observables: %d" % self.ny)

        dsetY = self.f.create_dataset(name="/measurements/y",
                                      shape=(self.num_conditions,
                                             self.num_timepoints, self.ny),
                                      chunks=(1, self.num_timepoints, self.ny),
                                      fillvalue=np.nan, dtype='f8',
                                      compression=self.compression)

        dsetSigmaY = self.f.create_dataset(name="/measurements/ysigma",
                                           shape=(self.num_conditions,
                                                  self.num_timepoints,
                                                  self.ny),
                                           chunks=(
                                               1, self.num_timepoints,
                                               self.ny),
                                           fillvalue=np.nan, dtype='f8',
                                           compression=self.compression)

        self.writeMeasurements(dsetY, dsetSigmaY)
        self.f.flush()

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
        for p in self.sbml_model.getListOfParameters():
            p = p.getId()
            if p.startswith('observable_') and not p.endswith('_sigma'):
                # observableName = p[len('observable_'):]
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
        # TODO: at the moment we cannot deal with replicate measurements
        #  (same condition, same observable, same timepoint) since not
        #  supported by AMICI

        # create inverse mapping for faster lookup
        condition_id_to_index = {
            name: idx for idx, name in enumerate(self.condition_ids)}
        observable_id_to_index = {
            name: idx for idx, name in enumerate(self.observable_ids)}

        # write measurements from each row in measurementDf
        for index, row in self.measurement_df.iterrows():
            condition_idx = condition_id_to_index[row.simulationConditionId]
            observable_idx = observable_id_to_index[row.observableId]
            time_idx = self.timepoints.index(row.time)

            if not np.isnan(dsetY[condition_idx, time_idx, observable_idx]):
                raise RuntimeError('Multiple measurements provided for '
                                   f'{row.simulationConditionId} - {row.observableId} time {row.time}'
                                   'Replicate measurements cannot be handled currently.')

            dsetY[condition_idx, time_idx, observable_idx] = \
                float(row['measurement'])

        # set sigmas
        # numeric sigma overrides have been mapped by petab already
        # will be the same for all timepoints
        # TODO: what if not? checked by petab?

        model_parameter_ids = self.model.getParameterIds()

        for condition_idx, condition_map in enumerate(self.parameter_mapping):
            for model_parameter_idx, mapped_parameter in enumerate(
                    condition_map):
                if isinstance(mapped_parameter, str):
                    continue

                model_parameter_id = model_parameter_ids[model_parameter_idx]
                if not model_parameter_id.startswith('noiseParameter1'):
                    continue

                observable_id = model_parameter_id[len('noiseParameter1_'):]
                observable_idx = observable_id_to_index[observable_id]
                dsetSigmaY[condition_idx, :, observable_idx] = mapped_parameter

    def generateHierarchicalOptimizationData(self):
        """
        Deal with offsets, proportionality factors and sigmas for hierarchical optimization

        Generate the respective index lists and mappings
        """

        # TODO: should check at the end of the function if scalingIndices lists
        # are non-overlapping

        observables = petab.get_observables(self.sbml_model, remove=False)
        sigmas = petab.get_sigmas(self.sbml_model, remove=False)
        print(observables, sigmas)

        observable_parameter_override_id_to_placeholder_id = {}
        noise_parameter_override_id_to_placeholder_id = {}

        for _, row in self.measurement_df.iterrows():
            # we trust that the number of overrides matches
            overrides = petab.split_parameter_replacement_list(
                row.observableParameters)
            placeholders = petab.get_placeholders(
                observables['observable_' + row.observableId]['formula'],
                row.observableId, 'observable')
            assert (len(overrides) == len(placeholders))
            for override, placeholder in zip(overrides, placeholders):
                if override in observable_parameter_override_id_to_placeholder_id:
                    if observable_parameter_override_id_to_placeholder_id[
                        override] != placeholder:
                        print(
                            f"WARNING: {override} overrides both {placeholder} and {
                            observable_parameter_override_id_to_placeholder_id[
                                override]}. Doubting that this will work...")
                        # TODO: append
                observable_parameter_override_id_to_placeholder_id[
                    override] = placeholder

            overrides = petab.split_parameter_replacement_list(
                row.noiseParameters)
            placeholders = petab.get_placeholders(
                sigmas['observable_' + row.observableId], row.observableId,
                'noise')
            assert (len(overrides) == len(placeholders))
            for override, placeholder in zip(overrides, placeholders):
                if override in noise_parameter_override_id_to_placeholder_id:
                    if noise_parameter_override_id_to_placeholder_id[
                        override] != placeholder:
                        print(
                            f"WARNING: {override} overrides both {placeholder} and {
                            noise_parameter_override_id_to_placeholder_id[
                                override]}. Doubting that this will work...")
                        # TODO: append
                noise_parameter_override_id_to_placeholder_id[
                    override] = placeholder

        print(observable_parameter_override_id_to_placeholder_id)
        print(noise_parameter_override_id_to_placeholder_id)

        offset_candidates = []
        scaling_candidates = []
        sigma_candidates = []

        for optimization_parameter_id in self.parameter_df.index[
            self.parameter_df.hierarchicalOptimization == 1]:
            # check which model parameter this one overrides

            # check in which observables this parameter occurs
            if optimization_parameter_id in observable_parameter_override_id_to_placeholder_id:
                placeholder_id = \
                observable_parameter_override_id_to_placeholder_id[
                    optimization_parameter_id]
                observable_id = '_'.join(placeholder_id.split('_')[1:])
                observable_formula = \
                observables['observable_' + observable_id]['formula']

                """
                print('optimization_parameter_id', optimization_parameter_id)
                print('placeholder_id', placeholder_id)
                print('observable_id', observable_id)
                print('observable_formula', observable_formula)
                """

                if petab.parameter_is_offset_parameter(placeholder_id,
                                                       observable_formula):
                    offset_candidates.append(optimization_parameter_id)
                elif petab.parameter_is_scaling_parameter(placeholder_id,
                                                          observable_formula):
                    scaling_candidates.append(optimization_parameter_id)
                else:
                    raise RuntimeError(
                        f'Parameter {optimization_parameter_id} selected for hierarchical optimization but is neither offset, proportionality or sigma parameter. Dunno what to do.')
            elif optimization_parameter_id in noise_parameter_override_id_to_placeholder_id:
                # TODO: what is there to check? formula - sigma == 0!
                sigma_candidates.append(optimization_parameter_id)
            else:
                # TODO: should also allow parameters which are no overrides
                # TODO ensure this is only output parameter
                raise RuntimeError(
                    f'Parameter {optimization_parameter_id} selected for hierarchical optimization but is neither offset, proportionality or sigma parameter. Dunno what to do.')

        print("offset_candidates:", offset_candidates)
        print("scaling_candidates:", scaling_candidates)
        print("sigma_candidates:", sigma_candidates)

        self.handleProportionalityFactors(scaling_candidates)
        self.handleOffsetParameter(
            offset_candidates)  # must call after handleProportionalityFactors
        self.handleSigmas(sigma_candidates)

    def handleOffsetParameter(self, offset_candidates):
        """
        Write list of offset parameters selected for hierarchical optimization
        """

        # don't create dataset if it would be empty
        if not len(offset_candidates):
            return

        # find indices for names
        offsetsForHierarchicalIndices = [
            self.optimization_parameter_name_to_index[x] for x in
            offset_candidates]

        print("Number of offset parameters for hierarchical optimization: %d" %
              len(offsetsForHierarchicalIndices))

        self.f.require_dataset("/offsetParameterIndices",
                               shape=(len(offsetsForHierarchicalIndices),),
                               dtype='<i4',
                               data=offsetsForHierarchicalIndices)

        # find usages for the selected parameters
        use = self.getAnalyticalParameterTable(offset_candidates, 'observable')

        self.f.require_dataset("/offsetParametersMapToObservables",
                               shape=(len(use), 3),
                               dtype='<i4', data=use)

    def getAnalyticalParameterTable(self, parametersForHierarchical, type):
        """Generate (scalingIdx, conditionIdx, observableIdx) table for all occurrences of the given parameter names

        Returns:
        list of (scalingIdx, conditionIdx, observableIdx) tuples
        """
        use = []
        for index, row in self.measurement_df.iterrows():
            # TODO : must handle sigmas separately
            if type == 'observable':
                overrides = petab.split_parameter_replacement_list(
                    row.observableParameters)
            elif type == 'noise':
                overrides = petab.split_parameter_replacement_list(
                    row.noiseParameters)
            else:
                raise ValueError(
                    "type must be noise or observable, but got" + type)

            for s in overrides:
                # print(s, parametersForHierarchical)
                try:
                    scalingIdx = parametersForHierarchical.index(s)
                    conditionIdx = self.condition_ids.index(
                        row.simulationConditionId)
                    observableIdx = self.observable_ids.index(row.observableId)
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

    def handleProportionalityFactors(self, scaling_candidates):
        """
        Write datasets specifying which proportionality factors to consider for
        hierarchical optimization
        """

        if not len(scaling_candidates):
            return

        scalingsForHierarchicalIndices = [
            self.optimization_parameter_name_to_index[x] for x in
            scaling_candidates]

        self.f.require_dataset("/scalingParameterIndices",
                               shape=(len(scalingsForHierarchicalIndices),),
                               dtype='<i4',
                               data=scalingsForHierarchicalIndices)
        print(
            "Number of proportionality factors for hierarchical optimization: %d" % len(
                scalingsForHierarchicalIndices))

        # find usages for the selected parameters
        use = self.getAnalyticalParameterTable(scaling_candidates,
                                               'observable')

        self.f.require_dataset("/scalingParametersMapToObservables",
                               shape=(len(use), 3),
                               dtype='<i4', data=use)

    def handleSigmas(self, sigma_candidates):
        """
        Write data for dealing with sigma parameters in hierarchical optimization
        """

        if not len(sigma_candidates):
            return

        sigmasForHierarchicalIndices = [
            self.optimization_parameter_name_to_index[x] for x in
            sigma_candidates]

        self.f.require_dataset("/sigmaParameterIndices",
                               shape=(len(sigmasForHierarchicalIndices),),
                               dtype='<i4',
                               data=sigmasForHierarchicalIndices)
        print("Number of sigmas for hierarchical optimization: %d" %
              len(sigmasForHierarchicalIndices))

        # find usages for the selected parameters
        use = self.getAnalyticalParameterTable(sigma_candidates, 'noise')

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
            '/amiciOptions/sens_ind', shape=(numParameters,), dtype="<i4",
            data=range(numParameters))

        self.f.require_dataset(
            '/amiciOptions/ts', shape=(len(self.timepoints),), dtype="f8",
            data=self.timepoints)

        # set parameter scaling: all log10, except for offsets which can be negative
        # ... for AMICI model parameters:
        offsetIndices = [i for i, p in enumerate(
            self.f['/parameters/modelParameterNames']) if
                         p.startswith('offset_')]
        pscale = np.full(numParameters, 2)
        pscale[offsetIndices] = 0
        self.f.require_dataset('/amiciOptions/pscale',
                               shape=(numParameters,), dtype="<i4",
                               data=pscale)

        # ... for all parameters for hierarchical optimization
        offsetIndices = [i for i, p in enumerate(
            self.f['/parameters/parameterNames']) if p.startswith('offset_')]
        # self.getAnalyticallyComputedSimulationParameterIndices()
        linParametersAmiciIndices = offsetIndices
        numOptimizationParameters = self.f['/parameters/parameterNames'].shape[
            0]
        self.f.require_dataset('/parameters/pscale',
                               shape=(numOptimizationParameters,), dtype="<i4",
                               data=[2 * (ip not in linParametersAmiciIndices)
                                     for ip in
                                     range(numOptimizationParameters)])

    def getAnalyticallyComputedSimulationParameterIndices(self):
        """
        Get model parameter index (not optimization parameter index) of all analytically computed parameters
        """
        parameterNamesModel = []
        if '/offsetParameterIndices' in self.f:
            parameterNamesOptimization = self.f['/parameters/parameterNames'][
                self.f['/offsetParameterIndices']]
            parameterNamesModel.extend(
                set([self.getGenericParameterName(o) for o in
                     parameterNamesOptimization]))

        if '/scalingParameterIndices' in self.f:
            parameterNamesOptimization = self.f['/parameters/parameterNames'][
                self.f['/scalingParameterIndices']]
            parameterNamesModel.extend(
                set([self.getGenericParameterName(o) for o in
                     parameterNamesOptimization]))

        if '/sigmaParameterIndices' in self.f:
            parameterNamesOptimization = self.f['/parameters/parameterNames'][
                self.f['/sigmaParameterIndices']]
            parameterNamesModel.extend(
                set([self.getGenericParameterName(o) for o in
                     parameterNamesOptimization]))

        return [self.f['/parameters/modelParameterNames'][:].tolist().index(p)
                for p in set(parameterNamesModel)]

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
            '/optimizationOptions/randomStarts',
            [numParams, numStartingPoints], 'f8')
        lower = self.f['/parameters/lowerBound'][:]
        upper = self.f['/parameters/upperBound'][:]
        startingPoints[:] = np.transpose(
            np.random.rand(numStartingPoints, numParams) * (
                    upper - lower) + lower)


def parse_cli_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        description='Create AMICI model and data for steadystate example.')

    parser.add_argument('-s', '--sbml', dest='sbml_file_name',
                        required=True,
                        help='SBML model filename (PEtab format)')

    parser.add_argument('-d', '--model-dir', dest='model_dir',
                        help='Model directory containing the python module')

    parser.add_argument('-n', '--model-name', dest='model_name',
                        required=True,
                        help='Name of the AMICI model module')

    parser.add_argument('-m', '--measurements', dest='measurement_file_name',
                        required=True,
                        help='Name of measurement table (PEtab format)')

    parser.add_argument('-c', '--conditions', dest='condition_file_name',
                        required=True,
                        help='Condition table (PEtab format)')

    parser.add_argument('-p', '--parameters', dest='parameter_file_name',
                        required=True,
                        help='Condition table (PEtab format)')

    parser.add_argument('-o', dest='hdf5_file_name', default='data.h5',
                        help='Name of HDF5 file to generate')

    # TODO read bounds and so on from there
    # parser.add_argument('-p', dest='parameter_file_name',
    #                    help='Parameter table', required=True)

    args = parser.parse_args()

    return args


def main():
    args = parse_cli_args()

    h5gen = HDF5DataGenerator(
        args.sbml_file_name,
        args.measurement_file_name,
        args.condition_file_name,
        args.parameter_file_name,
        args.model_dir, args.model_name)
    h5gen.generateFile(args.hdf5_file_name)


if __name__ == "__main__":
    main()
