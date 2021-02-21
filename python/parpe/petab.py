"""Functions related to handling of PEtab files"""
from numbers import Number
from typing import Tuple, Dict, List

from petab import C as ptc
import petab
import pandas as pd


def get_parameter_override_id_to_placeholder_id(
        observable_df: pd.DataFrame,
        measurement_df: pd.DataFrame
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """Get map of parameter overrides from PEtab measurement table to
    placeholder parameter IDs from observable table that they are overriding.

    Arguments:
        observable_df: PEtab observable table
        measurement_df: PEtab measurement table

    Returns:
        Dictionaries mapping parameter overrides to placeholder parameter IDs
        that they are overriding, for observable parameters and noise
        parameters.
    """

    observable_placeholders = \
        get_observable_placeholders_by_observable(observable_df)
    noise_placeholders = \
        get_noise_placeholders_by_observable(observable_df)

    observable_parameter_override_id_to_placeholder_id = {}
    noise_parameter_override_id_to_placeholder_id = {}

    def _add(override_to_placeholders):
        """Add current overrides to destination map"""
        assert (len(overrides) == len(placeholders))

        for override, placeholder in zip(overrides, placeholders):
            if isinstance(override, Number):
                # TODO: we ignore numeric overrides, cannot currently use for
                #  hierarchical optimization
                continue

            try:
                override_to_placeholders[override].append(placeholder)
            except KeyError:
                override_to_placeholders[override] = [placeholder]

    for obs_id, obs_par, noise_par \
            in zip(measurement_df[ptc.OBSERVABLE_ID],
                   measurement_df[ptc.OBSERVABLE_PARAMETERS],
                   measurement_df[ptc.NOISE_PARAMETERS]):
        # observable parameters
        overrides = petab.split_parameter_replacement_list(obs_par)
        placeholders = observable_placeholders[obs_id]
        _add(observable_parameter_override_id_to_placeholder_id)
        # noise parameters
        overrides = petab.split_parameter_replacement_list(noise_par)
        placeholders = noise_placeholders[obs_id]
        assert (len(overrides) == len(placeholders))
        _add(noise_parameter_override_id_to_placeholder_id)

    return observable_parameter_override_id_to_placeholder_id, \
           noise_parameter_override_id_to_placeholder_id


def get_observable_placeholders_by_observable(
        observable_df: pd.DataFrame) -> Dict[str, List[str]]:
    """Get observable placeholder parameters for all observables

    Arguments:
        observable_df: PEtab observable table

    Returns:
        Dictionary mapping observable ID to list of placeholder parameters
    """

    return {obs_id: petab.get_formula_placeholders(formula, obs_id, 'observable')
            for obs_id, formula
            in zip(observable_df.index, observable_df[ptc.OBSERVABLE_FORMULA])}


def get_noise_placeholders_by_observable(
        observable_df: pd.DataFrame) -> Dict[str, List[str]]:
    """Get noise placeholder parameters for all observables

    Arguments:
        observable_df: PEtab observable table

    Returns:
        Dictionary mapping observable ID to list of placeholder parameters
    """

    return {obs_id: petab.get_formula_placeholders(formula, obs_id, 'noise')
            for obs_id, formula
            in zip(observable_df.index, observable_df[ptc.NOISE_FORMULA])}
