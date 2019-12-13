#!/usr/bin/env python3
"""
Generate HDF5 file for parPE with fixed parameters and measurements for an
AMICI-imported SBML model based on tables with fixed parameters and training
data in the PEtab format https://github.com/ICB-DCM/PEtab
"""

import numpy as np
import h5py
import sys
import argparse
import petab
import importlib
import amici
import pandas as pd
import libsbml
from pandas import DataFrame
from colorama import init as init_colorama
from colorama import Fore
from numbers import Number
from typing import Iterable, Any, Collection, Optional, Dict, List, Tuple


def to_float_if_float(x: Any) -> Any:
    try:
        return float(x)
    except ValueError:
        return x


def unique_ordered(seq: Iterable) -> list:
    """
    Make unique, preserving order of first occurrence

    Arguments:
        seq: any sequence
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def requires_preequilibration(measurement_df: DataFrame) -> bool:
    return 'preequilibrationConditionId' in measurement_df \
            and not np.issubdtype(
                measurement_df.preequilibrationConditionId.dtype, np.number)


def write_string_array(f: h5py.File,
                       path: str,
                       strings: Collection):
    """
    Write string array to hdf5

    Arguments:
        f: h5py.File
        path: path of the dataset to create
        strings: list of strings
    """
    dt = h5py.special_dtype(vlen=str)
    dset = f.create_dataset(path, (len(strings),), dtype=dt)
    dset[:] = [s.encode('utf8') for s in strings]
    f.flush()


def write_float_array(f, path, values, dtype='f8'):
    """
    Write float array to hdf5

    Arguments:
        f: h5py.File
        path: path of the dataset to create
        values: array to write
        dtype: datatype
    """
    dset = f.create_dataset(path, (len(values),), dtype=dtype)
    dset[:] = values
    f.flush()


def write_int_array(f, path, values, dtype='<i4'):
    """
    Write integer array to hdf5

    Arguments:
        f: h5py.File
        path: path of the dataset to create
        values: array to write
        dtype: datatype
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
        measurement_df: PEtab measurement table
        condition_df: PEtab condition table
        parameter_df: PEtab parameter table
        sbml_model: SBML model
        condition_ids: numpy.array condition IDs (different condition vectors,
            both simulation and preequilibration)
        num_condition_vectors:
            Number of condition vectors, including simulation and
            preequilibration. Not necessarily equal to the number of
            simulations.
        unique_timepoints: time points for which there is data
        f: h5py.File
            hdf5 file which is being created
    """

    def __init__(self, sbml_file, measurement_file, condition_file,
                 parameter_file, model_output_dir, model_name, verbose=1):
        """
        fileNameSBML: filename of model SBML file (PEtab-style)
        fileMeasurements: filename of measurement file
        fileFixedParameters: filename with AMICI fixed parameter vectors for
            all conditions referred to in the measurement file
        fileParameters: PEtab parameter table filename
        """
        self.condition_ids = None
        self.f: Optional[h5py.File] = None
        self.num_condition_vectors: int = 0
        self.unique_timepoints = None
        self.parameter_mapping: Optional[petab.ParMappingDict] = None
        self.parameter_scale_mapping: Optional[petab.ScaleMappingDict] = None
        self.optimization_parameter_name_to_index: Dict[str, int] = {}
        self.nk: int = 0
        self.condition_map = None
        self.observable_ids = None
        self.ny: int = 0

        self.verbose = verbose

        self.petab_problem = petab.Problem.from_files(
            sbml_file=sbml_file,
            condition_file=condition_file,
            measurement_file=measurement_file,
            parameter_file=parameter_file
        )
        # ensure we have valid inputs
        petab.lint_problem(self.petab_problem)

        self.sbml_model: libsbml.Model = self.petab_problem.sbml_model
        self.condition_df: pd.DataFrame = self.petab_problem.condition_df
        self.measurement_df: pd.DataFrame = self.petab_problem.measurement_df
        self.parameter_df: pd.DataFrame = self.petab_problem.parameter_df

        # parse data
        self.parse_measurement_file()
        self.parse_fixed_parameters_file()

        # load model
        sys.path.insert(0, model_output_dir)
        model_module = importlib.import_module(model_name)
        self.amici_model = model_module.getModel()

        # index for no reference/preequilibration condition
        self.NO_PREEQ_CONDITION_IDX: int = -1

        # value for unmapped model parameter in opt<->sim mapping
        self.UNMAPPED_PARAMETER: int = -1

        # scriptDir = os.path.dirname(os.path.realpath(__file__))

        # hdf5 dataset compression
        self.compression = "gzip"

    def parse_measurement_file(self) -> None:
        """
        Read cost_fun file and determine number of conditions and timepoints
        """

        if self.verbose:
            print(Fore.CYAN + "Measurements shape:", self.measurement_df.shape)

        # Get list of used conditions
        self.condition_ids = unique_ordered(
            self.measurement_df.simulationConditionId.values)
        if requires_preequilibration(self.measurement_df):
            self.condition_ids = unique_ordered(
                [*self.condition_ids,
                 *self.measurement_df.preequilibrationConditionId.values])
        self.num_condition_vectors = len(self.condition_ids)

        # when using adjoint sensitivities, we cannot keep inf
        # -> consider late timepoint as steady-state
        print(Fore.GREEN + "Changing t = Inf to t = 1e8, since we cannot use "
                           "Inf with adjoints.")
        self.measurement_df.loc[
            self.measurement_df.time == np.inf, 'time'] = 1e8

        # list of unique timepoints, just for info and AMICI model default
        # setting
        self.unique_timepoints = sorted(self.measurement_df.time.unique())

        print(Fore.CYAN + "Num condition vectors: ",
              self.num_condition_vectors)
        print(Fore.CYAN + "Num timepoints: ", self.unique_timepoints,
              len(self.unique_timepoints))

    def parse_fixed_parameters_file(self) -> None:
        """
        Load file and select relevant conditions
        """
        print(Fore.CYAN + "Fixed parameters orginal: ",
              self.condition_df.shape)

        # drop conditions that do not have measurements
        drop_rows = [label for label in self.condition_df.index
                     if label not in self.condition_ids]
        self.condition_df.drop(drop_rows, 0, inplace=True)

        print(Fore.CYAN + "Fixed parameters usable: ",
              self.condition_df.shape)

    def generate_file(self, hdf5_file_name) -> None:
        """
        Create the output file

        fileNameH5:   filename of HDF5 file that is to be generated
        """
        self.f = h5py.File(hdf5_file_name, "w")

        self.save_metadata()

        print(Fore.GREEN + "Generating simulation condition list...")
        self.generate_simulation_condition_map()

        print(Fore.GREEN + "Generating parameter list...")
        self.generate_parameter_list()

        print(Fore.GREEN + "Generating fixed parameters matrix...")
        self.generate_fixed_parameter_matrix()

        print(Fore.GREEN + "Generating measurement matrix...")
        self.generate_measurement_matrices()

        print(Fore.GREEN + "Handling scaling parameters...")
        self.generate_hierarchical_optimization_data()

        print(Fore.GREEN + "Copying default AMICI options...")
        self.copy_amici_options()

        print(Fore.GREEN + "Writing default optimization options...")
        self.write_optimization_options()

    def save_metadata(self):
        """Save some extra information in the generated file"""

        g = self.f.require_group('/metadata')

        g.attrs['invocation'] = ' '.join(sys.argv)
        g.attrs['amici_version'] = amici.__version__
        # TODO: parPE version

    def generate_parameter_list(self) -> None:
        """
        Optimization to simulation parameter mapping. Write parameter names.
        """

        # simulation parameters from model
        model_parameter_ids = np.array(self.amici_model.getParameterIds())
        write_string_array(self.f, "/parameters/modelParameterNames",
                           model_parameter_ids)
        print(Fore.CYAN + "Number of model parameters:",
              len(model_parameter_ids))

        print(Fore.CYAN + "Number of optimization parameters:",
              len(self.parameter_df))
        write_string_array(self.f, "/parameters/parameterNames",
                           self.parameter_df.index.values[
                               (self.parameter_df.estimate == 1)
                           & ~self.parameter_df.index.isin(
                        self.amici_model.getFixedParameterIds())])

        self.generate_simulation_to_optimization_parameter_mapping()

        self.f.flush()

    def generate_simulation_to_optimization_parameter_mapping(self) -> None:
        """
        Create dataset n_parameters_simulation x n_conditions with indices of
        respective parameters in pararameters_optimization
        """

        # get list of tuple of parameters dicts for condition
        self.parameter_mapping = self.petab_problem \
            .get_optimization_to_simulation_parameter_mapping(
                warn_unmapped=False)
        self.parameter_scale_mapping = \
            petab.get_optimization_to_simulation_scale_mapping(
                mapping_par_opt_to_par_sim=self.parameter_mapping,
                parameter_df=self.petab_problem.parameter_df,
                measurement_df=self.petab_problem.measurement_df
            )

        print(self.parameter_mapping)

        amici_model_parameter_ids = self.amici_model.getParameterIds()

        # Merge and preeq and sim parameters, filter fixed parameters
        for condition_idx, \
            ((condition_map_preeq, condition_map_sim),
             (condition_scale_map_preeq, condition_scale_map_sim)) \
                in enumerate(zip(self.parameter_mapping,
                                 self.parameter_scale_mapping)):
            petab.merge_preeq_and_sim_pars_condition(
                condition_map_preeq, condition_map_sim,
                condition_scale_map_preeq, condition_scale_map_sim,
                condition_idx)

            # PEtab mapping may contain fixed parameters, filter out first
            # NOTE: parameter ordering might differ from AMICI ordering
            self.parameter_mapping[condition_idx] = np.array(
                [condition_map_sim[par]
                 for par in amici_model_parameter_ids])
            self.parameter_scale_mapping[condition_idx] = np.array(
                [condition_scale_map_sim[par]
                 for par in amici_model_parameter_ids])

        # Translate parameter ID mapping to index mapping
        # create inverse mapping for faster lookup
        optimization_parameter_name_to_index = {
            name: idx for idx, name
            in enumerate(
                self.parameter_df.index[
                    (self.parameter_df.estimate == 1)
                    & (~self.parameter_df.index.isin(
                        self.amici_model.getFixedParameterIds()))
                    ])}
        # print(optimization_parameter_name_to_index)
        self.optimization_parameter_name_to_index = \
            optimization_parameter_name_to_index

        # use in-memory matrix, don't write every entry to file directly
        num_model_parameters = self.amici_model.np()
        mapping_matrix = np.zeros(
            shape=(num_model_parameters, self.condition_map.shape[0]),
            dtype='<i4')
        override_matrix = np.full(shape=mapping_matrix.shape,
                                  fill_value=np.nan)
        # for each condition index vector
        for condition_idx, par_map \
                in enumerate(self.parameter_mapping):

            # for each model parameter
            for model_parameter_idx, mapped_parameter \
                    in enumerate(par_map):
                mapped_parameter = to_float_if_float(mapped_parameter)
                try:
                    (mapped_idx, override) = self.get_index_mapping_for_par(
                        mapped_parameter, optimization_parameter_name_to_index)
                    mapping_matrix[model_parameter_idx, condition_idx] = \
                        mapped_idx
                    override_matrix[model_parameter_idx, condition_idx] = \
                        override
                except IndexError as e:
                    print(Fore.RED + "Error in parameter mapping:", e)
                    print(model_parameter_idx, mapped_parameter)
                    print(self.parameter_mapping)
                    raise e

        # write to file
        write_parameter_map(self.f, mapping_matrix, override_matrix,
                            num_model_parameters, self.compression)
        write_scale_map(self.f, self.parameter_scale_mapping,
                        self.parameter_df, self.amici_model)

    def get_index_mapping_for_par(
            self, mapped_parameter: Any,
            optimization_parameter_name_to_index: Dict[str, int]
    ) -> Tuple[int, float]:
        if isinstance(mapped_parameter, str):
            # actually a mapped optimization parameter
            try:
                return (optimization_parameter_name_to_index[mapped_parameter],
                        np.nan)
            except KeyError:
                # This is a fixed parameter which is to be replaced
                # by nominalValue
                return (self.UNMAPPED_PARAMETER,
                        self.parameter_df.loc[mapped_parameter,
                                              'nominalValue'])
        elif np.isnan(mapped_parameter):
            # This condition does not use any parameter override.
            # We override with 0.0, NAN will cause AMICI warnings.
            return self.UNMAPPED_PARAMETER, 1.0
        else:
            # Have numeric override
            return self.UNMAPPED_PARAMETER, mapped_parameter

    def generate_fixed_parameter_matrix(self) -> None:
        """
        Write fixed parameters dataset (nFixedParameters x nConditions).
        """

        fixed_parameter_ids = self.amici_model.getFixedParameterIds()
        self.nk = len(fixed_parameter_ids)
        print(Fore.CYAN + "Number of fixed parameters:",
              len(fixed_parameter_ids))

        # Create in-memory table, write all at once for speed
        fixed_parameter_matrix = np.full(
            shape=(self.nk, self.num_condition_vectors),
            fill_value=np.nan)
        for i in range(len(fixed_parameter_ids)):
            self.handle_fixed_parameter(i, fixed_parameter_ids[i],
                                        fixed_parameter_matrix)

        self.create_fixed_parameter_dataset_and_write_attributes(
            fixed_parameter_ids, fixed_parameter_matrix)

        self.f.flush()

    def create_fixed_parameter_dataset_and_write_attributes(
            self, fixed_parameter_ids: Collection[str], data) -> h5py.Dataset:
        """
        Create fixed parameters data set and annotations
        """
        self.f.require_group("/fixedParameters")

        write_string_array(self.f, "/fixedParameters/parameterNames",
                           fixed_parameter_ids)
        write_string_array(self.f, "/fixedParameters/conditionNames",
                           self.condition_ids)

        # chunked for reading experiment-wise
        nk = len(fixed_parameter_ids)

        if len(data):
            dset = self.f.create_dataset("/fixedParameters/k",
                                         dtype='f8', chunks=(nk, 1),
                                         compression=self.compression,
                                         data=data)
            # set dimension scales
            dset.dims.create_scale(
                self.f['/fixedParameters/parameterNames'], 'parameterNames')
            dset.dims.create_scale(
                self.f['/fixedParameters/conditionNames'], 'conditionNames')
            dset.dims[0].attach_scale(
                self.f['/fixedParameters/parameterNames'])
            dset.dims[1].attach_scale(
                self.f['/fixedParameters/conditionNames'])
        else:
            dset = self.f.create_dataset("/fixedParameters/k",
                                         dtype='f8')

        return dset

    def generate_simulation_condition_map(self):
        """
        Write index map for independent simulations to be performed
        (preequilibrationConditionIdx, simulationConditionIdx) referencing the
        fixed parameter table.

        Parse preequilibrationConditionId column, and write to hdf5 dataset.
        If conditionRef is empty, set to -1, otherwise to condition index.

        """
        simulations = \
            self.petab_problem.get_simulation_conditions_from_measurement_df()
        condition_map = np.full(shape=(len(simulations), 2),
                                fill_value=self.NO_PREEQ_CONDITION_IDX)
        condition_id_to_idx = {cid: idx for idx, cid
                               in enumerate(self.condition_ids)}
        if simulations.shape[1] == 2:
            # preeq always nan, we only need simulation condition id
            condition_map[:, 1] = [condition_id_to_idx[condition_id]
                                   for condition_id in simulations.iloc[:, 0]]
        else:
            # we need to set both preeq and simulation condition
            # (mind different ordering preeq/sim sim/preeq here and in PEtab
            for sim_idx, (sim_id, preeq_id) \
                    in enumerate(simulations.iloc[:, 0:2].values):
                condition_map[sim_idx] = [condition_id_to_idx[preeq_id],
                                          condition_id_to_idx[sim_id]]

        print(Fore.CYAN + "Number of simulation conditions:",
              len(simulations))

        self.condition_map = condition_map
        self.f.create_dataset("/fixedParameters/simulationConditions",
                              dtype="<i4",
                              data=condition_map)

    def handle_fixed_parameter(self, parameter_index: int, parameter_name: str,
                               dset: h5py.Dataset) -> None:
        """
        Extract parameter values from data table or model and set to dset

        Lookup order is condition table > model parameter > model species

        condition table is expected to have "globalized names", i.e.
        reactionId_localParameterId

        NOTE: parameters which are not provided in fixedParametersDf or model
        are set to 0.0

        Arguments:
            parameter_index:
                Index of the current fixed parameter
            parameter_name:
                Name of the current parameter
            dset:
                2D array-like (nFixedParameters x nConditions) where the
                parameter value is to be set
        """

        if parameter_name in self.condition_df.columns:
            # Parameter in condition table
            dset[parameter_index, :] = \
                self.condition_df.loc[self.condition_ids,
                                      parameter_name].values
            return

        sbml_parameter = self.sbml_model.getParameter(parameter_name)
        if sbml_parameter:
            # Parameter value from model
            dset[parameter_index, :] = sbml_parameter.getValue()
            return

        # Is this a species?
        sbml_species = self.sbml_model.getSpecies(parameter_name)
        if sbml_species:
            # A constant species might have been turned in to a model
            # parameter
            # TODO: we dont do any conversion here, although we would
            #  want to have concentration currently there is only 1.0
            dset[parameter_index, :] = \
                sbml_species.getInitialConcentration() \
                if sbml_species.isSetInitialConcentration() \
                else sbml_species.getInitialAmount()
            return

        # We need to check for "globalized" parameter names too
        # (reactionId_localParameterId)
        # model has localParameterId, data file has globalized name
        global_name = get_global_name_for_local_parameter(
            self.sbml_model, parameter_name)
        if global_name:
            sbml_parameter = self.sbml_model.getParameter(
                global_name)
            if sbml_parameter:
                # Use model parameter value
                dset[parameter_index, :] = \
                    sbml_parameter.getValue()
                return

        print(Fore.YELLOW
              + "Warning: Fixed parameter not "
                "found in ExpTable, setting to 0.0: ",
              parameter_name)
        dset[parameter_index, :] = 0.0

    def generate_measurement_matrices(self):
        """
        Generate matrices with training data for conditions, observables and
        timepoints. Matrices are numbered by simulation conditions, each of
        dimension of the respective num_timepoints x num_observables.
        num_observables will always be ny from the model, since this is what is
        required for AMICI. num_timepoints is condition dependent
        """
        if petab.measurement_table_has_timepoint_specific_mappings(
                self.measurement_df):
            raise RuntimeError("Timepoint-specific overrides are not yet "
                               "supported.")

        self.f.create_group("/measurements")
        self.observable_ids = self.amici_model.getObservableIds()
        # trim observable_ TODO: should be done in amici import
        self.observable_ids = [o[len('observable_'):]
                               if o.startswith('observable_') else o
                               for o in self.observable_ids]
        self.ny = self.amici_model.ny
        write_string_array(self.f, "/measurements/observableNames",
                           self.observable_ids)

        print(Fore.CYAN + "Number of observables:", self.ny)

        self.write_measurements()
        self.f.flush()

    def write_measurements(self):
        """
        Write measurements to hdf5 dataset
        """

        # create inverse mapping for faster lookup
        observable_id_to_index = {
            name: idx for idx, name in enumerate(self.observable_ids)}

        for sim_idx, (preeq_cond_idx, sim_cond_idx) \
                in enumerate(self.condition_map):
            # print("Condition", sim_idx, (preeq_cond_idx, sim_cond_idx))

            row_filter = 1
            row_filter &= self.measurement_df.simulationConditionId \
                == self.condition_ids[sim_cond_idx]
            if preeq_cond_idx == self.NO_PREEQ_CONDITION_IDX:
                row_filter &= \
                    self.measurement_df.preequilibrationConditionId.isnull()
            else:
                row_filter &= \
                    self.measurement_df.preequilibrationConditionId \
                    == self.condition_ids[preeq_cond_idx]
            cur_mes_df = self.measurement_df.loc[row_filter, :]
            if not len(cur_mes_df):
                raise AssertionError("No measurements for this condition: ",
                                     sim_idx, (preeq_cond_idx, sim_cond_idx))

            # get the required, possibly non-unique (replicates) timepoints
            grp = cur_mes_df.groupby(['observableId', 'time']).size() \
                .reset_index().rename(columns={0:'sum'}) \
                .groupby(['time'])['sum'].agg(['max']).reset_index()
            timepoints = np.sort(np.repeat(grp['time'], grp['max'])).tolist()
            mes = np.full(shape=(len(timepoints), self.ny),
                           fill_value=np.nan)
            sd = mes.copy()

            # write measurements from each row in measurementDf
            for index, row in cur_mes_df.iterrows():
                observable_idx = observable_id_to_index[row.observableId]
                time_idx = timepoints.index(row.time)
                while not np.isnan(mes[time_idx, observable_idx]):
                    time_idx += 1
                    if timepoints[time_idx] != row.time:
                        raise AssertionError(
                            "Replicate handling failed for "
                            f'{row.simulationConditionId} - {row.observableId}'
                            f' time {row.time}\n' + str(cur_mes_df))
                mes[time_idx, observable_idx] = float(row['measurement'])
                sigma = to_float_if_float(row['noiseParameters'])
                if isinstance(sigma, float):
                    sd[time_idx, observable_idx] = sigma

            self.f.create_dataset(name=f"/measurements/y/{sim_idx}",
                                  data=mes, dtype='f8',
                                  compression=self.compression)
            self.f.create_dataset(name=f"/measurements/ysigma/{sim_idx}",
                                  data=sd, dtype='f8',
                                  compression=self.compression)
            self.f.create_dataset(name=f"/measurements/t/{sim_idx}",
                                  data=timepoints, dtype='f8',
                                  compression=self.compression)

    def get_parameter_override_id_to_placeholder_id(self, observables, sigmas):
        # override -> [observable|sigma ids]
        observable_parameter_override_id_to_placeholder_id = {}
        noise_parameter_override_id_to_placeholder_id = {}

        # TODO: this seems to be redundant with the creation of the mapping
        #  tables later on
        for _, row in self.measurement_df.iterrows():
            # we trust that the number of overrides matches

            # observable parameters
            overrides = petab.split_parameter_replacement_list(
                row.observableParameters)
            placeholders = petab.get_placeholders(
                observables['observable_' + row.observableId]['formula'],
                row.observableId, 'observable')
            assert (len(overrides) == len(placeholders))
            for override, placeholder in zip(overrides, placeholders):
                if isinstance(override, Number):
                    # ignore numeric overrides, cannot currently use for
                    # hierarchical optimization TODO
                    continue

                try:
                    observable_parameter_override_id_to_placeholder_id[
                        override].append(placeholder)
                except KeyError:
                    observable_parameter_override_id_to_placeholder_id[
                        override] = [placeholder]

            # noise parameters
            overrides = petab.split_parameter_replacement_list(
                row.noiseParameters)
            placeholders = petab.get_placeholders(
                sigmas['observable_' + row.observableId], row.observableId,
                'noise')
            assert (len(overrides) == len(placeholders))
            for override, placeholder in zip(overrides, placeholders):
                if isinstance(override, Number):
                    continue

                try:
                    noise_parameter_override_id_to_placeholder_id[
                        override].append(placeholder)
                except KeyError:
                    noise_parameter_override_id_to_placeholder_id[
                        override] = [placeholder]

        return observable_parameter_override_id_to_placeholder_id, \
               noise_parameter_override_id_to_placeholder_id

    def generate_hierarchical_optimization_data(self, verbose=1):
        """
        Deal with offsets, proportionality factors and sigmas for hierarchical
        optimization

        Generate the respective index lists and mappings
        """

        if 'hierarchicalOptimization' not in self.parameter_df:
            print(Fore.YELLOW + 'Missing hierarchicalOptimization column in '
                  'parameter table. Skipping.')
            return

        observables = petab.get_observables(self.sbml_model, remove=False)
        sigmas = petab.get_sigmas(self.sbml_model, remove=False)

        if verbose:
            print(Fore.CYAN + "Observables:")
            print(Fore.CYAN, observables)
            print()
            print(Fore.CYAN + "Sigmas:")
            print(Fore.CYAN, sigmas)
            print()

        observable_parameter_override_id_to_placeholder_id, \
        noise_parameter_override_id_to_placeholder_id = \
            self.get_parameter_override_id_to_placeholder_id(observables,
                                                             sigmas)

        # print(Fore.CYAN, observable_parameter_override_id_to_placeholder_id)
        # print(Fore.CYAN, noise_parameter_override_id_to_placeholder_id)

        offset_candidates = set()
        scaling_candidates = set()
        sigma_candidates = set()

        hierarchical_candidates = self.parameter_df.index[
            (self.parameter_df.estimate == 1)
            & (self.parameter_df.hierarchicalOptimization == 1)]

        for optimization_parameter_id in hierarchical_candidates:
            # check which model parameter this one overrides

            # check in which observables this parameter occurs
            if optimization_parameter_id \
                    in observable_parameter_override_id_to_placeholder_id:
                placeholder_ids = \
                    observable_parameter_override_id_to_placeholder_id[
                        optimization_parameter_id]
                for placeholder_id in placeholder_ids:
                    observable_id = '_'.join(placeholder_id.split('_')[1:])
                    observable_formula = \
                        observables['observable_' + observable_id]['formula']

                    """
                    print('optimization_parameter_id', optimization_parameter_id)
                    print('placeholder_id', placeholder_id)
                    print('observable_id', observable_id)
                    print('observable_formula', observable_formula)
                    """

                    if petab.parameter_is_offset_parameter(
                            placeholder_id, observable_formula):
                        offset_candidates.add(optimization_parameter_id)
                    elif petab.parameter_is_scaling_parameter(
                            placeholder_id, observable_formula):
                        scaling_candidates.add(optimization_parameter_id)
                    else:
                        raise RuntimeError(
                            f'Parameter {optimization_parameter_id} selected '
                            'for hierarchical optimization but is neither '
                            'offset, proportionality or sigma parameter. '
                            f'Dunno what to do. {observable_id}')
            elif optimization_parameter_id \
                    in noise_parameter_override_id_to_placeholder_id:
                # TODO: what is there to check? formula - sigma == 0!
                sigma_candidates.add(optimization_parameter_id)
            else:
                # TODO: should also allow parameters which are no overrides
                # TODO ensure this is only output parameter
                raise RuntimeError(
                    f'Parameter {optimization_parameter_id} selected '
                    'for hierarchical optimization but is neither '
                    'offset, proportionality or sigma parameter. '
                    'Dunno what to do.')

        print(Fore.CYAN + "offset_candidates:", offset_candidates)
        print(Fore.CYAN + "scaling_candidates:", scaling_candidates)
        print(Fore.CYAN + "sigma_candidates:", sigma_candidates)

        # check if scalingIndices lists are non-overlapping
        for x in offset_candidates:
            if x in scaling_candidates:
                raise RuntimeError(
                    f"Determined {x} as candidate for both offset and scaling.")
            if x in sigma_candidates:
                raise RuntimeError(
                    f"Determined {x} as candidate for both offset and sigma.")
        for x in scaling_candidates:
            if x in sigma_candidates:
                raise RuntimeError(
                    f"Determined {x} as candidate for both scaling and sigma.")

        offset_candidates = list(offset_candidates)
        scaling_candidates = list(scaling_candidates)
        sigma_candidates = list(sigma_candidates)

        self.handle_proportionality_factors(scaling_candidates)
        # must call after handle_proportionality_factors
        self.handle_offset_parameter(offset_candidates)
        self.handle_sigmas(sigma_candidates)

        # check:
        if '/scalingParametersMapToObservables' in self.f \
                and '/sigmaParametersMapToObservables' in self.f:
            scaling_map = self.f['/scalingParametersMapToObservables'][:]
            sigma_map = self.f['/sigmaParametersMapToObservables'][:]
            df = pd.DataFrame(data=dict(scaling_id=scaling_map[:, 0],
                                        condition_id=scaling_map[:, 1],
                                        observable_id=scaling_map[:, 2]))
            df.set_index(['observable_id', 'condition_id'], inplace=True)

            df2 = pd.DataFrame(data=dict(sigma_id=sigma_map[:, 0],
                                         condition_id=sigma_map[:, 1],
                                         observable_id=sigma_map[:, 2]))
            df2.set_index(['observable_id', 'condition_id'], inplace=True)
            df = df.join(df2)
            del df2

            # TODO: smarter check
            if df.isnull().values.any():
                print(Fore.YELLOW + "Couldn't verify that parameter selection "
                                    "for hierarchical optimization is ok.")
            else:
                df_grouped = \
                    df.groupby(['scaling_id', 'sigma_id']).size().reset_index()
                # must be the same, otherwise one scaling is used with
                # multiple sigma
                if len(df_grouped) != len(df_grouped.scaling_id.unique()):
                    raise AssertionError("Scaling parameter selected for hierarchical "
                                         "optimization is used with multiple sigmas.")
                # TODO: same check for offsets

    def handle_offset_parameter(self, offset_candidates):
        """
        Write list of offset parameters selected for hierarchical optimization
        """

        # don't create dataset if it would be empty
        if not len(offset_candidates):
            return

        # find indices for names
        offsets_for_hierarchical_indices = [
            self.optimization_parameter_name_to_index[x] for x in
            offset_candidates]
        order = np.argsort(offsets_for_hierarchical_indices)
        offsets_for_hierarchical_indices = \
            [offsets_for_hierarchical_indices[i] for i in order]
        offset_candidates = [offset_candidates[i] for i in order]

        print(Fore.CYAN + "Number of offset parameters for hierarchical "
                          "optimization: %d"
              % len(offsets_for_hierarchical_indices))

        self.f.require_dataset("/offsetParameterIndices",
                               shape=(len(offsets_for_hierarchical_indices),),
                               dtype='<i4',
                               data=offsets_for_hierarchical_indices)

        # find usages for the selected parameters
        use = self.get_analytical_parameter_table(offset_candidates,
                                                  'observable')

        self.f.require_dataset("/offsetParametersMapToObservables",
                               shape=(len(use), 3),
                               dtype='<i4', data=use)

    def get_analytical_parameter_table(
            self,
            hierarchical_candidate_ids: list,
            parameter_type: str) -> list:
        """Generate (scalingIdx, conditionIdx, observableIdx) table for all
        occurrences of the given parameter names.

        Parameters:
            hierarchical_candidate_ids: Ids of optimization parameters for
                hierarchical optimization. This table depends on ordering of
                this list.
            parameter_type:
                'observable' or 'noise'

        Returns:
        list of (scalingIdx, conditionIdx, observableIdx) tuples
        """

        condition_id_to_index = {name: idx for idx, name in
                                 enumerate(self.condition_ids)}
        # need list, not ndarray
        condition_map_list = [list(x) for x in self.condition_map]

        use = []
        for index, row in self.measurement_df.iterrows():
            # TODO : must handle sigmas separately
            if parameter_type == 'observable':
                overrides = petab.split_parameter_replacement_list(
                    row.observableParameters)
            elif parameter_type == 'noise':
                overrides = petab.split_parameter_replacement_list(
                    row.noiseParameters)
            else:
                raise ValueError(
                    "type must be noise or observable, but got" + parameter_type)

            sim_cond_idx = \
                condition_id_to_index[row.simulationConditionId]
            preeq_cond_idx = self.NO_PREEQ_CONDITION_IDX
            if not np.isnan(row.preequilibrationConditionId):
                preeq_cond_idx = condition_id_to_index[
                    row.preequilibrationConditionId]

            for s in overrides:
                #print(s, parametersForHierarchical)
                try:
                    scalingIdx = hierarchical_candidate_ids.index(s)
                except ValueError:
                    continue # current parameter not in list

                conditionIdx = condition_map_list.index(
                    [preeq_cond_idx, sim_cond_idx])
                observableIdx = self.observable_ids.index(row.observableId)
                tup = (scalingIdx, conditionIdx, observableIdx)

                # Don't add a new line for each timepoint
                # We don't allow separate parameters for individual time-points
                # (Can be implemented via different observables)
                if not tup in use:
                    use.append(tup)

        if not len(use):
            raise AssertionError("Candidates were: " + str(hierarchical_candidate_ids) + " but nothing usable found")

        return use

    def handle_proportionality_factors(self, scaling_candidates):
        """
        Write datasets specifying which proportionality factors to consider for
        hierarchical optimization
        """

        if not len(scaling_candidates):
            return

        scalingsForHierarchicalIndices = [
            self.optimization_parameter_name_to_index[x] for x in
            scaling_candidates]
        order = np.argsort(scalingsForHierarchicalIndices)
        scalingsForHierarchicalIndices = \
            [scalingsForHierarchicalIndices[i] for i in order]
        scaling_candidates = [scaling_candidates[i] for i in order]


        self.f.require_dataset("/scalingParameterIndices",
                               shape=(len(scalingsForHierarchicalIndices),),
                               dtype='<i4',
                               data=scalingsForHierarchicalIndices)
        print(Fore.CYAN, "Number of proportionality factors for "
                         "hierarchical optimization: %d"
              % len(scalingsForHierarchicalIndices))

        # find usages for the selected parameters
        use = self.get_analytical_parameter_table(scaling_candidates,
                                                  'observable')

        self.f.require_dataset("/scalingParametersMapToObservables",
                               shape=(len(use), 3),
                               dtype='<i4', data=use)

    def handle_sigmas(self, sigma_candidates):
        """
        Write data for dealing with sigma parameters in hierarchical optimization

        Parameters:
            sigma_candidates:
                IDs of optimization parameters which are sigmas to be
                hierarchically computed. No particular order.
        """

        if not len(sigma_candidates):
            return

        # must sort by indices for C++ code
        sigmas_for_hierarchical_indices = [
            self.optimization_parameter_name_to_index[x] for x in
            sigma_candidates]
        order = np.argsort(sigmas_for_hierarchical_indices)
        sigmas_for_hierarchical_indices = \
            [sigmas_for_hierarchical_indices[i] for i in order]
        sigma_candidates = [sigma_candidates[i] for i in order]

        self.f.require_dataset("/sigmaParameterIndices",
                               shape=(len(sigmas_for_hierarchical_indices),),
                               dtype='<i4',
                               data=sigmas_for_hierarchical_indices)
        print(Fore.CYAN + "Number of sigmas for hierarchical optimization: %d"
              % len(sigmas_for_hierarchical_indices))

        # find usages for the selected parameters
        use = self.get_analytical_parameter_table(sigma_candidates, 'noise')

        self.f.require_dataset("/sigmaParametersMapToObservables",
                               shape=(len(use), 3),
                               dtype='<i4', data=use)

    def copy_amici_options(self):
        """
        Write simulation options
        """
        g = self.f.require_group("/amiciOptions")
        g.attrs['sensi'] = amici.SensitivityOrder_first
        g.attrs['sensi_meth'] = amici.SensitivityMethod_adjoint
        g.attrs['tstart'] = 0.0
        g.attrs['atol'] = 1e-14
        g.attrs['interpType'] = amici.InterpolationType_hermite
        g.attrs['ism'] = amici.InternalSensitivityMethod_simultaneous
        g.attrs['iter'] = amici.NonlinearSolverIteration_newton
        g.attrs['linsol'] = amici.LinearSolver_KLU
        g.attrs['lmm'] = amici.LinearMultistepMethod_BDF
        g.attrs['maxsteps'] = 10000
        # fail fast to retry with looser tolerances
        g.attrs['maxstepsB'] = 5000
        g.attrs['newton_preeq'] = 0
        g.attrs['nmaxevent'] = 0
        g.attrs['ordering'] = 0
        g.attrs['rtol'] = 1e-6
        g.attrs['stldet'] = 1

        num_model_parameters = \
            self.f['/parameters/modelParameterNames'].shape[0]

        # parameter indices w.r.t. which to compute sensitivities
        self.f.require_dataset(
            '/amiciOptions/sens_ind', shape=(num_model_parameters,),
            dtype="<i4",
            data=range(num_model_parameters))

        # TODO not really meaningful anymore - remove?
        self.f.require_dataset(
            '/amiciOptions/ts', shape=(len(self.unique_timepoints),), dtype="f8",
            data=self.unique_timepoints)

    def get_analytically_computed_optimization_parameter_indices(self):
        """
        Get optimization parameter index of all analytically computed
        parameters.
        """
        indices = []
        if '/offsetParameterIndices' in self.f:
            indices.extend(self.f['/offsetParameterIndices'])

        if '/scalingParameterIndices' in self.f:
            indices.extend(self.f['/scalingParameterIndices'])

        if '/sigmaParameterIndices' in self.f:
            indices.extend(self.f['/sigmaParameterIndices'])

        return list(set(indices))

    def write_optimization_options(self):
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

        self.write_bounds()
        self.write_starting_points()

    def write_bounds(self):
        """
        Parameter bounds for optimizer

        Offset parameters are allowed to be negative
        """
        optimized_par_df = \
            self.parameter_df.loc[self.parameter_df.estimate == 1
                                  & (~self.parameter_df.index.isin(
                        self.amici_model.getFixedParameterIds())), :]
        self.f.require_dataset('/parameters/lowerBound',
                               shape=optimized_par_df.lowerBound.shape,
                               data=optimized_par_df.lowerBound, dtype='f8')
        self.f.require_dataset('/parameters/upperBound',
                               shape=optimized_par_df.upperBound.shape,
                               data=optimized_par_df.upperBound, dtype='f8')

    def write_starting_points(self):
        """
        Write a list of random starting points sampled as specified in PEtab
        file.
        """
        num_params = self.f['/parameters/parameterNames'].shape[0]
        num_starting_points = 100
        np.random.seed(0)

        starting_points = self.f.require_dataset(
            '/optimizationOptions/randomStarts',
            [num_params, num_starting_points], 'f8')

        starting_points[:] = \
            self.petab_problem.sample_parameter_startpoints(
                num_starting_points)

        # Write nominal values for testing purposes
        if 'nominalValue' in self.parameter_df:
            self.f['/parameters/nominalValues'] = \
                self.parameter_df.nominalValue[
                    self.parameter_df.estimate == 1]


def petab_scale_to_amici_scale(scale_str):
    """Convert PEtab parameter scaling string to AMICI scaling integer"""

    if scale_str == 'lin':
        return amici.ParameterScaling_none
    if scale_str == 'log':
        return amici.ParameterScaling_ln
    if scale_str == 'log10':
        return amici.ParameterScaling_log10
    raise ValueError("Invalid pscale " + scale_str)


def get_global_name_for_local_parameter(sbml_model, needle_parameter_id):
    for reaction in sbml_model.getListOfReactions():
        kl = reaction.getKineticLaw()
        for p in kl.getListOfParameters():
            parameter_id = p.getId()
            if parameter_id.endswith(needle_parameter_id):
                return f'{reaction.getId()}_{parameter_id}'
    return None


def write_parameter_map(f: h5py.File, mapping_matrix: np.array,
                        override_matrix: np.array, num_model_parameters: int,
                        compression=None):
    """Write parameter mapping"""
    f.require_dataset(
        name='/parameters/parameterOverrides',
        chunks=(num_model_parameters, 1),
        dtype='f8',
        fillvalue=np.nan,
        compression=compression,
        shape=override_matrix.shape,
        data=override_matrix)

    f.require_dataset(
        name='/parameters/optimizationSimulationMapping',
        chunks=(num_model_parameters, 1),
        dtype='<i4',
        fillvalue=-1,
        compression=compression,
        shape=mapping_matrix.shape,
        data=mapping_matrix)


def write_scale_map(f: h5py.File, parameter_scale_mapping: List[List[str]],
                    parameter_df: pd.DataFrame, amici_model: amici.Model):
    """Write parameter scale mapping to HDF5 dataset"""

    # for simulation
    # set parameter scaling for all parameters
    pscale = np.zeros(shape=(len(parameter_scale_mapping),
                             len(parameter_scale_mapping[0])))
    for i, cond_scale_list in enumerate(parameter_scale_mapping):
        for j, s in enumerate(cond_scale_list):
            pscale[i, j] = petab_scale_to_amici_scale(s)

    f.require_dataset('/parameters/pscaleSimulation',
                      shape=pscale.shape,
                      dtype="<i4",
                      data=pscale)

    # for cost function parameters
    pscale = np.array([petab_scale_to_amici_scale(s)
                       for s in parameter_df.parameterScale.values[
                           (parameter_df.estimate == 1)
                           & ~parameter_df.index.isin(
                               amici_model.getFixedParameterIds())]])
    f.require_dataset('/parameters/pscaleOptimization',
                      shape=pscale.shape,
                      dtype="<i4",
                      data=pscale)


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

    args = parser.parse_args()

    return args


def main():
    init_colorama(autoreset=True)

    args = parse_cli_args()

    h5gen = HDF5DataGenerator(
        args.sbml_file_name,
        args.measurement_file_name,
        args.condition_file_name,
        args.parameter_file_name,
        args.model_dir, args.model_name)
    h5gen.generate_file(args.hdf5_file_name)


if __name__ == "__main__":
    main()
