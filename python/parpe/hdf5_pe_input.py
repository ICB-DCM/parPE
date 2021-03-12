"""Functions for generating parPE parameter estimation HDF5 input files"""
import argparse
import logging
import sys
from typing import Any, Collection, Optional, Dict, Tuple

import amici
import coloredlogs
import h5py
import numpy as np
import pandas as pd
import petab
from amici.petab_import import PREEQ_INDICATOR_ID
from amici.petab_import import petab_scale_to_amici_scale
from amici.petab_objective import subset_dict
from colorama import Fore
from colorama import init as init_colorama
from pandas import DataFrame
# noinspection PyPep8Naming
from petab import C as ptc
from petab import to_float_if_float

from .hdf5 import write_string_array
from .hierarchical_optimization import (get_candidates_for_hierarchical,
                                        get_analytical_parameter_table)
from .misc import get_amici_model, unique_ordered

logger = logging.getLogger(__name__)

# TODO: transpose all datasets?


# For condition mapping table in HDF5 file: value indicating no
# reference/preequilibration condition
NO_PREEQ_CONDITION_IDX: int = -1

# For parameter mapping table in HDF5 file: value indicating unmapped model
# parameter in opt<->sim mapping
UNMAPPED_PARAMETER: int = -1


def requires_preequilibration(measurement_df: DataFrame) -> bool:
    return ptc.PREEQUILIBRATION_CONDITION_ID in measurement_df \
           and not np.issubdtype(
        measurement_df[ptc.PREEQUILIBRATION_CONDITION_ID].dtype,
        np.number)


class HDF5DataGenerator:
    """
    Generate HDF5 file with fixed parameters and measurements for an
    AMICI-imported SBML model based on a PEtab problem.

    Attributes:
        amici_model: AMICI model for which to estimate parameters
        petab_problem: PEtab optimization problem (will be modified in place)
        compression: h5py compression to be used
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

    def __init__(self,
                 petab_problem: petab.Problem,
                 amici_model: amici.Model,
                 verbose=1):
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
        self.petab_problem: petab.Problem = petab_problem
        self.amici_model: amici.Model = amici_model

        # ensure we have valid inputs
        petab.lint_problem(self.petab_problem)

        # index for no reference/preequilibration condition
        self.NO_PREEQ_CONDITION_IDX: int = NO_PREEQ_CONDITION_IDX

        # value for unmapped model parameter in opt<->sim mapping
        self.UNMAPPED_PARAMETER: int = UNMAPPED_PARAMETER

        # hdf5 dataset compression
        self.compression = "gzip"

    def generate_file(self, hdf5_file_name: str) -> None:
        """
        Create the output file

        Arguments:
            hdf5_file_name: filename of HDF5 file that is to be generated
        """
        self.f = h5py.File(hdf5_file_name, "w")

        self._save_metadata()
        self._collect_simulation_conditions()
        self._process_condition_table()

        logger.info(Fore.GREEN + "Generating simulation condition list...")
        self._generate_simulation_condition_map()

        logger.info(Fore.GREEN + "Generating parameter list...")
        self._generate_parameter_list()
        self._generate_simulation_to_optimization_parameter_mapping()

        logger.info(Fore.GREEN + "Generating measurement matrix...")
        self._generate_measurement_matrices()

        logger.info(Fore.GREEN + "Handling scaling parameters...")
        self._generate_hierarchical_optimization_data()

        logger.info(Fore.GREEN + "Copying default AMICI options...")
        self._write_amici_options()

        logger.info(Fore.GREEN + "Writing default optimization options...")
        write_optimization_options(self.f)
        self._write_bounds()
        self._write_starting_points()

    def _collect_simulation_conditions(self) -> None:
        """
        Read cost_fun file and determine number of conditions and timepoints
        """
        measurement_df = self.petab_problem.measurement_df
        if self.verbose:
            logger.info(f"{Fore.CYAN}Number of data points:{Fore.RESET} "
                        f"{measurement_df.shape[0]}")

        # Get list of used conditions vectors
        self.condition_ids = unique_ordered(
            measurement_df[ptc.SIMULATION_CONDITION_ID].values)

        if requires_preequilibration(measurement_df):
            self.condition_ids = unique_ordered(
                [*self.condition_ids,
                 *measurement_df[
                     ptc.PREEQUILIBRATION_CONDITION_ID].values])
        self.num_condition_vectors = len(self.condition_ids)

        self.condition_id_to_index = {name: idx for idx, name in
                                      enumerate(self.condition_ids)}

        # list of unique timepoints, just for info and AMICI model default
        # setting
        self.unique_timepoints = sorted(measurement_df[ptc.TIME].unique())

        logger.info(f"{Fore.CYAN}Num condition vectors:{Fore.RESET} "
                    f"{self.num_condition_vectors}")
        logger.info(f"{Fore.CYAN}Num timepoints:{Fore.RESET} "
                    f"{len(self.unique_timepoints)}")
        logger.debug(f"Timepoints: {self.unique_timepoints}")

    def _process_condition_table(self) -> None:
        """
        Load file and select relevant conditions
        """
        condition_df = self.petab_problem.condition_df
        if ptc.PARAMETER_NAME in condition_df:
            # we don't need that, so we drop it that we don't need to care
            # later on
            condition_df.drop(columns=ptc.PARAMETER_NAME)

        logger.info(f"{Fore.CYAN}Conditions original:{Fore.RESET} "
                    f"{condition_df.shape}")

        # drop conditions that do not have measurements
        drop_rows = [label for label in condition_df.index
                     if label not in self.condition_ids]
        condition_df.drop(index=drop_rows, inplace=True)
        logger.info(f"{Fore.CYAN}Conditions with data:{Fore.RESET} "
                    f"{condition_df.shape}")

    def _save_metadata(self):
        """Save some extra information in the generated file"""

        g = self.f.require_group('/metadata')

        g.attrs['invocation'] = ' '.join(sys.argv)
        g.attrs['amici_version'] = amici.__version__
        g.attrs['petab_version'] = petab.__version__
        # TODO: parPE version
        # g.attrs['parpe_version'] = parpe.__version__

        # Model info
        # Allows for checking in C++ code whether we are likely to use the
        # correct model
        g = self.f.require_group('/model')
        g.attrs['model_name'] = self.amici_model.getName()
        write_string_array(g, "observableIds",
                           self.amici_model.getObservableIds())
        write_string_array(g, "parameterIds",
                           self.amici_model.getParameterIds())
        write_string_array(g, "fixedParameterIds",
                           self.amici_model.getFixedParameterIds())
        write_string_array(g, "stateIds",
                           self.amici_model.getStateIds())

    def _generate_parameter_list(self) -> None:
        """
        Optimization to simulation parameter mapping. Write parameter names.
        """

        # simulation parameters from model
        model_parameter_ids = np.array(self.amici_model.getParameterIds())
        # TODO: rename to "ID"
        write_string_array(self.f, "/parameters/modelParameterNames",
                           model_parameter_ids)
        logger.info(f"{Fore.CYAN}Number of model parameters: {Fore.RESET}"
                    f"{len(model_parameter_ids)}")

        self.problem_parameter_ids = self.petab_problem \
            .get_optimization_parameters()

        # sanity check: estimated parameters should not be AMICI fixed
        # parameters
        fixed_opt_pars = set(self.problem_parameter_ids) \
                         & set(self.amici_model.getFixedParameterIds())
        if fixed_opt_pars:
            raise RuntimeError(f"Parameter {fixed_opt_pars} are to be "
                               "optimized, but are fixed parameters in the "
                               "model. This should not happen.")

        logger.info(f"{Fore.CYAN}Number of optimization parameters: "
                    f"{Fore.RESET}{len(self.problem_parameter_ids)}")

        write_string_array(self.f, "/parameters/parameterNames",
                           self.problem_parameter_ids)

    def _generate_simulation_to_optimization_parameter_mapping(self) -> None:
        """
        Create dataset n_parameters_simulation x n_conditions with indices of
        respective parameters in parameters_optimization
        """

        # get list of tuple of parameters dicts for all conditions
        self.parameter_mapping = self.petab_problem \
            .get_optimization_to_simulation_parameter_mapping(
            warn_unmapped=False, scaled_parameters=False)

        variable_par_ids = self.amici_model.getParameterIds()
        fixed_par_ids = self.amici_model.getFixedParameterIds()

        # Translate parameter ID mapping to index mapping
        # create inverse mapping for faster lookup
        logger.debug(f"Problem parameters: {self.problem_parameter_ids}")
        optimization_parameter_name_to_index = {
            name: idx for idx, name in enumerate(self.problem_parameter_ids)}
        self.optimization_parameter_name_to_index = \
            optimization_parameter_name_to_index

        # use in-memory matrix, don't write every entry to file directly
        num_model_parameters = self.amici_model.np()
        mapping_matrix = np.zeros(
            shape=(num_model_parameters, self.condition_map.shape[0]),
            dtype='<i4')
        override_matrix = np.full(shape=mapping_matrix.shape,
                                  fill_value=np.nan)
        pscale_matrix = np.full(shape=mapping_matrix.shape, dtype='<i4',
                                fill_value=amici.ParameterScaling_none)

        # AMICI fixed parameters
        self.nk = len(fixed_par_ids)
        logger.info(f"{Fore.CYAN}Number of fixed parameters:{Fore.RESET} "
                    f"{len(fixed_par_ids)}")
        # TODO: can fixed parameter differ between same condition vector used
        #  in different sim x preeq combinations?
        fixed_parameter_matrix = np.full(
            shape=(self.nk, self.num_condition_vectors),
            fill_value=np.nan)

        # For handling initial states
        species_in_condition_table = [
            col for col in self.petab_problem.condition_df
            if self.petab_problem.sbml_model.getSpecies(col) is not None]

        # for reinitialization
        state_id_to_idx = {
            id: i for i, id in enumerate(self.amici_model.getStateIds())}
        state_idxs_reinitialization_all = []

        # Merge and preeq and sim parameters, filter fixed parameters
        for condition_idx, \
            (condition_map_preeq, condition_map_sim,
             condition_scale_map_preeq, condition_scale_map_sim) \
                in enumerate(self.parameter_mapping):

            preeq_cond_idx, sim_cond_idx = self.condition_map[condition_idx]
            preeq_cond_id = self.condition_ids[preeq_cond_idx] \
                if preeq_cond_idx != NO_PREEQ_CONDITION_IDX else None
            sim_cond_id = self.condition_ids[sim_cond_idx]
            logger.debug(f"preeq_cond_idx: {preeq_cond_idx}, "
                         f"preeq_cond_id: {preeq_cond_id}, "
                         f"sim_cond_idx: {sim_cond_idx}, "
                         f"sim_cond_id: {sim_cond_id}")

            if len(condition_map_preeq) != len(condition_scale_map_preeq) \
                    or len(condition_map_sim) != len(condition_scale_map_sim):
                raise AssertionError(
                    "Number of parameters and number of parameter "
                    "scales do not match.")
            if len(condition_map_preeq) \
                    and len(condition_map_preeq) != len(condition_map_sim):
                # logger.debug(
                #     f"Preequilibration parameter map: {condition_map_preeq}")
                # logger.debug(f"Simulation parameter map: {condition_map_sim}")
                raise AssertionError(
                    "Number of parameters for preequilibration "
                    "and simulation do not match.")

            # state indices for reinitialization
            #  for current simulation condition
            state_idxs_for_reinitialization_cur = []

            # TODO: requires special handling of initial concentrations
            if species_in_condition_table:
                # set indicator fixed parameter for preeq
                # (we expect here, that this parameter was added during AMICI
                # model import and that it was not added by the user with a
                # different meaning...)
                if preeq_cond_idx != NO_PREEQ_CONDITION_IDX:
                    condition_map_preeq[PREEQ_INDICATOR_ID] = 1.0
                    condition_scale_map_preeq[PREEQ_INDICATOR_ID] = ptc.LIN

                condition_map_sim[PREEQ_INDICATOR_ID] = 0.0
                condition_scale_map_sim[PREEQ_INDICATOR_ID] = ptc.LIN

                def _set_initial_concentration(condition_id, species_id,
                                               init_par_id,
                                               par_map, scale_map):
                    value = petab.to_float_if_float(
                        self.petab_problem.condition_df.loc[
                            condition_id, species_id])
                    if isinstance(value, float):
                        # numeric initial state
                        par_map[init_par_id] = value
                        scale_map[init_par_id] = ptc.LIN
                    else:
                        # parametric initial state
                        try:
                            # try find in mapping
                            par_map[init_par_id] = par_map[value]
                            scale_map[init_par_id] = scale_map[value]
                        except KeyError:
                            # otherwise look up in parameter table
                            if (self.petab_problem.parameter_df.loc[
                                value, ptc.ESTIMATE] == 0):
                                par_map[init_par_id] = \
                                    self.petab_problem.parameter_df.loc[
                                        value, ptc.NOMINAL_VALUE]
                            else:
                                par_map[init_par_id] = value

                            if (ptc.PARAMETER_SCALE
                                    not in self.petab_problem.parameter_df
                                    or not self.petab_problem.parameter_df.loc[
                                        value, ptc.PARAMETER_SCALE]):
                                scale_map[init_par_id] = ptc.LIN
                            else:
                                scale_map[init_par_id] = \
                                    self.petab_problem.parameter_df.loc[
                                        value, ptc.PARAMETER_SCALE]

                for species_id in species_in_condition_table:
                    # for preequilibration
                    init_par_id = f'initial_{species_id}_preeq'

                    # need to set dummy value for preeq parameter anyways, as it
                    #  is expected below (set to 0, not nan, because will be
                    #  multiplied with indicator variable in initial assignment)
                    condition_map_sim[init_par_id] = 0.0
                    condition_scale_map_sim[init_par_id] = ptc.LIN

                    if preeq_cond_idx != NO_PREEQ_CONDITION_IDX:
                        value = petab.to_float_if_float(
                            self.petab_problem.condition_df.loc[
                                preeq_cond_id, species_id])
                        if isinstance(value, float):
                            condition_map_sim[init_par_id] = value
                        _set_initial_concentration(
                            preeq_cond_id, species_id, init_par_id,
                            condition_map_preeq,
                            condition_scale_map_preeq)
                        # enable state reinitialization
                        self.f['/fixedParameters/simulationConditions'][
                            condition_idx, 2] = 1

                    # for simulation
                    init_par_id = f'initial_{species_id}_sim'
                    _set_initial_concentration(
                        sim_cond_id, species_id, init_par_id,
                        condition_map_sim,
                        condition_scale_map_sim)

                    # mark for reinitialization
                    species_idx = state_id_to_idx[species_id]
                    state_idxs_for_reinitialization_cur.append(species_idx)
            state_idxs_reinitialization_all.append(state_idxs_for_reinitialization_cur)
            logger.debug(f"condition_map_preeq: {condition_map_preeq}, "
                         f"condition_map_sim: {condition_map_sim}")

            # split into fixed and variable parameters:
            condition_map_preeq_var = condition_scale_map_preeq_var = None
            condition_map_preeq_fix = None
            if condition_map_preeq:
                condition_map_preeq_var, condition_map_preeq_fix = \
                    subset_dict(condition_map_preeq, variable_par_ids,
                                fixed_par_ids)
                condition_scale_map_preeq_var, _ = \
                    subset_dict(condition_scale_map_preeq, variable_par_ids,
                                fixed_par_ids)

            condition_map_sim_var, condition_map_sim_fix = \
                subset_dict(condition_map_sim, variable_par_ids, fixed_par_ids)
            condition_scale_map_sim_var, condition_scale_map_sim_fix = \
                subset_dict(condition_scale_map_sim, variable_par_ids,
                            fixed_par_ids)

            if condition_map_preeq:
                # merge after having removed potentially fixed parameters
                # which may differ between preequilibration and simulation
                petab.merge_preeq_and_sim_pars_condition(
                    condition_map_preeq_var, condition_map_sim_var,
                    condition_scale_map_preeq_var, condition_scale_map_sim_var,
                    condition_idx)

            # mapping for each model parameter
            for model_parameter_idx, model_parameter_id \
                    in enumerate(variable_par_ids):
                mapped_parameter = condition_map_sim[model_parameter_id]
                mapped_parameter = to_float_if_float(mapped_parameter)
                try:
                    (mapped_idx, override) = self.get_index_mapping_for_par(
                        mapped_parameter, optimization_parameter_name_to_index)
                    mapping_matrix[model_parameter_idx, condition_idx] = \
                        mapped_idx
                    override_matrix[model_parameter_idx, condition_idx] = \
                        override
                except IndexError as e:
                    logger.critical(f"{Fore.RED}Error in parameter mapping: "
                                    f"{Fore.RESET} {e}")
                    logger.critical(f"{model_parameter_idx}, "
                                    f"{mapped_parameter}")
                    logger.critical(f"{self.parameter_mapping}")
                    raise e

            # Set parameter scales for simulation
            pscale_matrix[:, condition_idx] = list(
                petab_scale_to_amici_scale(condition_scale_map_sim[amici_id])
                for amici_id in variable_par_ids)

            fixed_parameter_matrix[:, self.condition_map[condition_idx, 1]] = \
                np.array([condition_map_sim_fix[par_id]
                          for par_id in fixed_par_ids])

            if condition_map_preeq:
                fixed_parameter_matrix[:,
                self.condition_map[condition_idx, 0]] = \
                    np.array([condition_map_preeq_fix[par_id]
                              for par_id in fixed_par_ids])

        # write to file
        self.create_fixed_parameter_dataset_and_write_attributes(
            fixed_par_ids, fixed_parameter_matrix)

        self.f.require_dataset('/parameters/pscaleSimulation',
                               shape=pscale_matrix.T.shape, dtype="<i4",
                               data=pscale_matrix.T)

        write_parameter_map(self.f, mapping_matrix, override_matrix,
                            num_model_parameters, self.compression)

        # for cost function parameters
        petab_opt_par_scale = \
            self.petab_problem.get_optimization_parameter_scales()
        pscale_opt_par = np.array(
            list(petab_scale_to_amici_scale(petab_opt_par_scale[par_id])
                 for par_id in self.problem_parameter_ids))

        self.f.require_dataset('/parameters/pscaleOptimization',
                               shape=pscale_opt_par.shape, dtype="<i4",
                               data=pscale_opt_par)

        # Ragged array of state indices for reinitialization
        data_type = h5py.vlen_dtype(np.dtype('int32'))
        dset = self.f.create_dataset(
            '/fixedParameters/reinitializationIndices',
            shape=(self.condition_map.shape[0],),
            dtype=data_type)
        for i, state_idxs in enumerate(state_idxs_reinitialization_all):
            dset[i] = state_idxs

        self.f.flush()

    def get_index_mapping_for_par(
            self, mapped_parameter: Any,
            optimization_parameter_name_to_index: Dict[str, int]
    ) -> Tuple[int, float]:
        """Get index mapping for a given parameter

        Arguments:
            mapped_parameter: value mapped to some model parameter
            optimization_parameter_name_to_index: Dictionary mapping
                optimization parameter IDs to their position in the
                optimization parameter vector
        Returns:
            Index of parameter to map to, and value for override matrix
        """
        if isinstance(mapped_parameter, str):
            # actually a mapped optimization parameter
            return (optimization_parameter_name_to_index[mapped_parameter],
                    np.nan)

        if np.isnan(mapped_parameter):
            # This condition does not use any parameter override.
            # We override with 0.0, NAN will cause AMICI warnings.
            return self.UNMAPPED_PARAMETER, 1.0

        # Have numeric override
        return self.UNMAPPED_PARAMETER, mapped_parameter

    def create_fixed_parameter_dataset_and_write_attributes(
            self, fixed_parameter_ids: Collection[str], data) -> h5py.Dataset:
        """
        Create fixed parameters data set and annotations
        """
        g = self.f.require_group("fixedParameters")

        write_string_array(g, "parameterNames",
                           fixed_parameter_ids)
        write_string_array(g, "conditionNames",
                           self.condition_ids)

        # chunked for reading experiment-wise
        nk = len(fixed_parameter_ids)

        if len(data):
            dset = g.create_dataset("k", dtype='f8', chunks=(nk, 1),
                                    compression=self.compression, data=data)
            # set dimension scales
            dset.dims.create_scale(g['parameterNames'], 'parameterNames')
            dset.dims.create_scale(g['conditionNames'], 'conditionNames')
            dset.dims[0].attach_scale(g['parameterNames'])
            dset.dims[1].attach_scale(g['conditionNames'])
        else:
            dset = g.create_dataset("k", dtype='f8')

        return dset

    def _generate_simulation_condition_map(self):
        """
        Write index map for independent simulations to be performed
        (preequilibrationConditionIdx, simulationConditionIdx) referencing the
        fixed parameter table.

        Get PEtab simulation conditions, and write to hdf5 dataset.
        If conditionRef is empty, set to NO_PREEQ_CONDITION_IDX, otherwise to
        respective condition index.
        """
        # PEtab condition list
        simulations = \
            self.petab_problem.get_simulation_conditions_from_measurement_df()
        # empty dataset
        condition_map = np.full(shape=(len(simulations), 2),
                                fill_value=self.NO_PREEQ_CONDITION_IDX)
        condition_id_to_idx = {condition_id: idx for idx, condition_id
                               in enumerate(self.condition_ids)}

        if simulations.shape[1] == 1:
            # preeq always nan, we only need simulation condition id
            condition_map[:, 1] = \
                list(condition_id_to_idx[condition_id]
                     for condition_id
                     in simulations[ptc.SIMULATION_CONDITION_ID])
        else:
            # We need to set both preeq and simulation condition
            # Mind different ordering preeq/sim sim/preeq here and in PEtab!
            for sim_idx, (sim_id, preeq_id) \
                    in enumerate(simulations.iloc[:, 0:2].values):
                if preeq_id:
                    condition_map[sim_idx, 0] = condition_id_to_idx[preeq_id]
                condition_map[sim_idx, 1] = condition_id_to_idx[sim_id]

        logger.info(f"{Fore.CYAN}Number of simulation conditions: {Fore.RESET}"
                    f" {len(simulations)}")

        self.condition_map = condition_map

        # append third column for state reinitialization
        # TODO : drop third column
        _condition_map = np.zeros((condition_map.shape[0],
                                   condition_map.shape[1] + 1), )
        _condition_map[:, :-1] = condition_map
        self.f.create_dataset("/fixedParameters/simulationConditions",
                              dtype="<i4",
                              data=_condition_map)

    def _generate_measurement_matrices(self):
        """
        Generate matrices with training data for conditions, observables and
        timepoints. Matrices are numbered by simulation conditions, each of
        dimension of the respective num_timepoints x num_observables.
        num_observables will always be ny from the model, since this is what is
        required for AMICI. num_timepoints is condition dependent
        """

        if petab.measurement_table_has_timepoint_specific_mappings(
                self.petab_problem.measurement_df):
            raise RuntimeError("Timepoint-specific overrides are not yet "
                               "supported.")

        self.f.create_group("/measurements")
        self.observable_ids = self.amici_model.getObservableIds()
        self.ny = self.amici_model.ny
        write_string_array(self.f, "/measurements/observableNames",
                           self.observable_ids)

        logger.info(f"{Fore.CYAN}Number of observables:{Fore.RESET} {self.ny}")

        self.write_measurements()
        self.f.flush()

    def write_measurements(self):
        """
        Write measurements to hdf5 dataset
        """

        # create inverse mapping for faster lookup
        observable_id_to_index = {
            name: idx for idx, name in enumerate(self.observable_ids)}

        measurement_df = self.petab_problem.measurement_df
        if ptc.NOISE_PARAMETERS not in measurement_df:
            measurement_df[ptc.NOISE_PARAMETERS] = np.nan

        for sim_idx, (preeq_cond_idx, sim_cond_idx) \
                in enumerate(self.condition_map):
            # print("Condition", sim_idx, (preeq_cond_idx, sim_cond_idx))

            cur_mes_df = self._get_measurements_for_condition(
                sim_idx, preeq_cond_idx, sim_cond_idx)

            # get the required, possibly non-unique (replicates) timepoints
            grp = cur_mes_df.groupby([ptc.OBSERVABLE_ID, ptc.TIME]).size() \
                .reset_index().rename(columns={0: 'sum'}) \
                .groupby([ptc.TIME])['sum'].agg(['max']).reset_index()

            timepoints = np.sort(np.repeat(grp[ptc.TIME], grp['max'])).tolist()
            mes = np.full(shape=(len(timepoints), self.ny), fill_value=np.nan)
            sd = mes.copy()

            # write measurements from each row in measurementDf
            for index, row in cur_mes_df.iterrows():
                observable_idx = observable_id_to_index[row[ptc.OBSERVABLE_ID]]
                time_idx = timepoints.index(row.time)
                while not np.isnan(mes[time_idx, observable_idx]):
                    time_idx += 1
                    if timepoints[time_idx] != row.time:
                        raise AssertionError(
                            "Replicate handling failed for "
                            f'{row[ptc.SIMULATION_CONDITION_ID]} - '
                            f'{row[ptc.OBSERVABLE_ID]}'
                            f' time {row[ptc.TIME]}\n' + str(cur_mes_df))
                mes[time_idx, observable_idx] = float(row[ptc.MEASUREMENT])
                sigma = to_float_if_float(row[ptc.NOISE_PARAMETERS])
                if isinstance(sigma, float):
                    sd[time_idx, observable_idx] = sigma

            # write to file
            g = self.f.require_group("measurements")
            g.create_dataset(name=f"y/{sim_idx}", data=mes, dtype='f8',
                             compression=self.compression)
            g.create_dataset(name=f"ysigma/{sim_idx}", data=sd, dtype='f8',
                             compression=self.compression)
            g.create_dataset(name=f"t/{sim_idx}", data=timepoints, dtype='f8',
                             compression=self.compression)

    def _get_measurements_for_condition(self, sim_idx,
                                        preeq_cond_idx, sim_cond_idx):
        """Get subset of measurement table for the given condition"""
        measurement_df = self.petab_problem.measurement_df

        row_filter = measurement_df[ptc.SIMULATION_CONDITION_ID] \
                     == self.condition_ids[sim_cond_idx]
        if ptc.PREEQUILIBRATION_CONDITION_ID in measurement_df:
            if preeq_cond_idx == self.NO_PREEQ_CONDITION_IDX:
                row_filter &= \
                    measurement_df[ptc.PREEQUILIBRATION_CONDITION_ID].isnull()
            else:
                row_filter &= \
                    measurement_df[ptc.PREEQUILIBRATION_CONDITION_ID] \
                    == self.condition_ids[preeq_cond_idx]
        cur_mes_df = measurement_df.loc[row_filter, :]
        if not len(cur_mes_df):
            # Should have been filtered out before
            raise AssertionError("No measurements for this condition: ",
                                 sim_idx, (preeq_cond_idx, sim_cond_idx))
        return cur_mes_df

    def _generate_hierarchical_optimization_data(self, verbose=1):
        """
        Deal with offsets, proportionality factors and sigmas for hierarchical
        optimization

        Generate the respective index lists and mappings
        """

        parameter_df = self.petab_problem.parameter_df
        if 'hierarchicalOptimization' not in parameter_df:
            logger.warning("Missing hierarchicalOptimization column in "
                           "parameter table. Skipping.")
            return

        if verbose:
            logger.info(f"{Fore.CYAN}Observables:")
            logger.info(f"{self.petab_problem.observable_df}")

        offset_candidates, scaling_candidates, sigma_candidates = \
            get_candidates_for_hierarchical(
                measurement_df=self.petab_problem.measurement_df,
                observable_df=self.petab_problem.observable_df,
                parameter_df=self.petab_problem.parameter_df)

        logger.info(f"{Fore.CYAN}offset_candidates:{Fore.RESET} "
                    f"{offset_candidates}")
        logger.info(f"{Fore.CYAN}scaling_candidates:{Fore.RESET} "
                    f"{scaling_candidates}")
        logger.info(f"{Fore.CYAN}sigma_candidates:{Fore.RESET} "
                    f"{sigma_candidates}")

        self._handle_proportionality_factors(scaling_candidates)
        # must call after _handle_proportionality_factors
        self.handle_offset_parameter(offset_candidates)
        self._handle_sigmas(sigma_candidates)
        self._check_hierarchical_optimization_tables()

    def _check_hierarchical_optimization_tables(self):
        """Try spotting potential error in hierarchical optimization tables"""

        if '/scalingParametersMapToObservables' not in self.f:
            return
        if '/sigmaParametersMapToObservables' not in self.f:
            return

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
            logger.warning("Couldn't verify that parameter selection "
                           "for hierarchical optimization is ok.")
        else:
            df_grouped = \
                df.groupby(['scaling_id', 'sigma_id']).size().reset_index()
            # must be the same, otherwise one scaling is used with
            # multiple sigma
            if len(df_grouped) != len(df_grouped.scaling_id.unique()):
                raise AssertionError(
                    "Scaling parameter selected for hierarchical "
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

        logger.info(f"{Fore.CYAN}Number of offset parameters for hierarchical "
                    f"optimization:{Fore.RESET} "
                    f"{len(offsets_for_hierarchical_indices)}")

        self.f.require_dataset("/offsetParameterIndices",
                               shape=(len(offsets_for_hierarchical_indices),),
                               dtype='<i4',
                               data=offsets_for_hierarchical_indices)

        # find usages for the selected parameters
        use = get_analytical_parameter_table(
            offset_candidates, 'observable', self.condition_id_to_index,
            self.petab_problem.measurement_df, self.observable_ids,
            self.condition_map, self.NO_PREEQ_CONDITION_IDX)

        self.f.require_dataset("/offsetParametersMapToObservables",
                               shape=(len(use), 3),
                               dtype='<i4', data=use)

    def _handle_proportionality_factors(self, scaling_candidates):
        """
        Write datasets specifying which proportionality factors to consider for
        hierarchical optimization
        """

        if not len(scaling_candidates):
            return

        scalings_for_hierarchical_indices = [
            self.optimization_parameter_name_to_index[x] for x in
            scaling_candidates]
        order = np.argsort(scalings_for_hierarchical_indices)
        scalings_for_hierarchical_indices = \
            [scalings_for_hierarchical_indices[i] for i in order]
        scaling_candidates = [scaling_candidates[i] for i in order]

        self.f.require_dataset("/scalingParameterIndices",
                               shape=(len(scalings_for_hierarchical_indices),),
                               dtype='<i4',
                               data=scalings_for_hierarchical_indices)
        logger.info(f"{Fore.CYAN}Number of proportionality factors for "
                    f"hierarchical optimization:{Fore.RESET} "
                    f"{len(scalings_for_hierarchical_indices)}")

        # find usages for the selected parameters
        use = get_analytical_parameter_table(
            scaling_candidates, 'observable', self.condition_id_to_index,
            self.petab_problem.measurement_df, self.observable_ids,
            self.condition_map, self.NO_PREEQ_CONDITION_IDX
        )

        self.f.require_dataset("/scalingParametersMapToObservables",
                               shape=(len(use), 3),
                               dtype='<i4', data=use)

    def _handle_sigmas(self, sigma_candidates):
        """
        Write data for dealing with sigma parameters in hierarchical
        optimization

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
        logger.info(f"{Fore.CYAN}Number of sigmas for hierarchical "
                    f"optimization:{Fore.RESET} "
                    f"{len(sigmas_for_hierarchical_indices)}")

        # find usages for the selected parameters
        use = get_analytical_parameter_table(
            sigma_candidates, 'noise', self.condition_id_to_index,
            self.petab_problem.measurement_df, self.observable_ids,
            self.condition_map, self.NO_PREEQ_CONDITION_IDX
        )

        self.f.require_dataset("/sigmaParametersMapToObservables",
                               shape=(len(use), 3),
                               dtype='<i4', data=use)

    def _write_amici_options(self) -> None:
        """
        Write simulation options
        """
        g = self.f.require_group("/amiciOptions")
        g.attrs['sensi'] = amici.SensitivityOrder.first
        g.attrs['sensi_meth'] = amici.SensitivityMethod.adjoint
        g.attrs['sensi_meth_preeq'] = amici.SensitivityMethod.adjoint
        # TODO PEtab support: get from file
        g.attrs['tstart'] = 0.0
        g.attrs['atol'] = 1e-14
        g.attrs['interpType'] = amici.InterpolationType_hermite
        g.attrs['ism'] = amici.InternalSensitivityMethod_simultaneous
        g.attrs['iter'] = amici.NonlinearSolverIteration_newton
        g.attrs['linsol'] = amici.LinearSolver_KLU
        g.attrs['lmm'] = amici.LinearMultistepMethod_BDF
        g.attrs['maxsteps'] = 10000
        # fail fast to retry with looser tolerances
        g.attrs['maxstepsB'] = 10000
        g.attrs['newton_preeq'] = 0
        g.attrs['nmaxevent'] = 0
        g.attrs['ordering'] = 0
        g.attrs['rtol'] = 1e-6
        g.attrs['stldet'] = 1

        # Required to handle parametric initial concentrations where newton
        # solver would fail with KLU error
        g.attrs['steadyStateSensitivityMode'] = \
            amici.SteadyStateSensitivityMode_simulationFSA

        num_model_parameters = \
            self.f['/parameters/modelParameterNames'].shape[0]

        # parameter indices w.r.t. which to compute sensitivities
        self.f.require_dataset(
            '/amiciOptions/sens_ind', shape=(num_model_parameters,),
            dtype="<i4",
            data=range(num_model_parameters))

        # TODO not really meaningful anymore - remove?
        self.f.require_dataset(
            '/amiciOptions/ts', shape=(len(self.unique_timepoints),),
            dtype="f8",
            data=self.unique_timepoints)

    def _get_analytically_computed_optimization_parameter_indices(self):
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

    def _write_bounds(self):
        """
        Parameter bounds for optimizer

        Offset parameters are allowed to be negative
        """
        lower_bound = np.array([petab.scale(
            self.petab_problem.parameter_df.loc[par_id, ptc.LOWER_BOUND],
            self.petab_problem.parameter_df.loc[par_id, ptc.PARAMETER_SCALE])
            for par_id in self.problem_parameter_ids])
        upper_bound = np.array([petab.scale(
            self.petab_problem.parameter_df.loc[par_id, ptc.UPPER_BOUND],
            self.petab_problem.parameter_df.loc[par_id, ptc.PARAMETER_SCALE])
            for par_id in self.problem_parameter_ids])

        self.f.require_dataset('/parameters/lowerBound',
                               shape=lower_bound.shape,
                               data=lower_bound, dtype='f8')
        self.f.require_dataset('/parameters/upperBound',
                               shape=upper_bound.shape,
                               data=upper_bound, dtype='f8')

    def _write_starting_points(self):
        """
        Write a list of random starting points sampled as specified in PEtab
        file.
        """
        num_params = len(self.problem_parameter_ids)
        num_starting_points = 100
        np.random.seed(0)

        starting_points = self.f.require_dataset(
            '/optimizationOptions/randomStarts',
            [num_params, num_starting_points], 'f8')

        starting_points[:] = \
            self.petab_problem.sample_parameter_startpoints(
                num_starting_points).T

        # Write nominal values for testing purposes
        if 'nominalValue' in self.petab_problem.parameter_df:
            # TODO: remove estimated=0, which should be part of mapping tables
            self.f['/parameters/nominalValues'] = list(petab.map_scale(
                self.petab_problem.parameter_df[ptc.NOMINAL_VALUE][
                    self.problem_parameter_ids],
                self.petab_problem.parameter_df[ptc.PARAMETER_SCALE][
                    self.problem_parameter_ids]))


def write_parameter_map(f: h5py.File, mapping_matrix: np.array,
                        override_matrix: np.array, num_model_parameters: int,
                        compression=None) -> None:
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


def write_optimization_options(f: h5py.File) -> None:
    """
    Create groups and write some default optimization settings
    """

    # set common options
    g = f.require_group('optimizationOptions')
    g.attrs['optimizer'] = 0  # IpOpt
    g.attrs['retryOptimization'] = 1
    g.attrs['hierarchicalOptimization'] = 1
    g.attrs['numStarts'] = 1

    # set IpOpt options
    g = f.require_group('optimizationOptions/ipopt')
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
    g = f.require_group('optimizationOptions/fmincon')
    g.attrs['MaxIter'] = 100
    g.attrs["TolX"] = 1e-8
    g.attrs["TolFun"] = 0
    g.attrs["MaxFunEvals"] = 1e7
    g.attrs["algorithm"] = np.string_("interior-point")
    g.attrs["GradObj"] = np.string_("on")
    g.attrs["display"] = np.string_("iter")

    # set CERES options
    g = f.require_group('optimizationOptions/ceres')
    g.attrs['max_num_iterations'] = 100

    # set toms611/SUMSL options
    g = f.require_group('optimizationOptions/toms611')
    g.attrs['mxfcal'] = 1e8

    # TODO mini-batch options


def parse_cli_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
        description='Create HDF5 input file for parameter estimation with '
                    'parPE.')

    parser.add_argument('-d', '--model-dir', dest='model_dir',
                        help='Model directory containing the python module')

    parser.add_argument('-n', '--model-name', dest='model_name',
                        required=True,
                        help='Name of the AMICI model module')

    parser.add_argument('-y', '--yaml', dest='petab_yaml',
                        required=True,
                        help='PEtab problem yaml file')

    parser.add_argument('-o', dest='hdf5_file_name', default='data.h5',
                        help='Name of HDF5 file to generate')

    args = parser.parse_args()

    return args


def main():
    """Generate HDF5 file based on CLI options"""
    init_colorama(autoreset=True)
    coloredlogs.DEFAULT_LEVEL_STYLES['debug'] = dict(color='white',
                                                     faint='true')
    coloredlogs.install(level='DEBUG', logger=logger)

    args = parse_cli_args()

    petab_problem = petab.Problem.from_yaml(args.petab_yaml)
    amici_model = get_amici_model(model_name=args.model_name,
                                  model_dir=args.model_dir)
    h5gen = HDF5DataGenerator(
        petab_problem=petab_problem,
        amici_model=amici_model)
    h5gen.generate_file(args.hdf5_file_name)
