"""Functions related to hierarchical optimization

https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz581/5538985
"""

from typing import Tuple, Dict, List
import numpy as np
import pandas as pd
import petab.C as ptc
import sympy as sp
from numbers import Number
from numpy import isnan
from petab import split_parameter_replacement_list

from .petab import get_parameter_override_id_to_placeholder_id


def parameter_is_scaling_parameter(parameter: str, formula: str) -> bool:
    """
    Check if is scaling parameter.

    Arguments:
        parameter: Some identifier.
        formula: Some sympy-compatible formula.

    Returns:
        ``True`` if parameter ``parameter`` is a scaling parameter in formula
         ``formula``.
    """

    sym_parameter = sp.sympify(parameter)
    sym_formula = sp.sympify(formula)

    return sym_parameter not in (sym_formula / sym_parameter).free_symbols


def parameter_is_offset_parameter(parameter: str, formula: str) -> bool:
    """
    Check if is offset parameter.

    Arguments:
        parameter: Some identifier.
        formula: Some sympy-compatible formula.

    Returns:
         ``True`` if parameter ``parameter`` is an offset parameter with
         positive sign in formula ``formula``.
    """

    sym_parameter = sp.sympify(parameter)
    sym_formula = sp.sympify(formula)

    return sym_parameter not in (sym_formula - sym_parameter).free_symbols


def get_candidates_for_hierarchical(
        observable_df: pd.DataFrame,
        measurement_df: pd.DataFrame,
        parameter_df: pd.DataFrame):
    """Based on PEtab files, check which parameters are suitable for
    hierarchical optimization.

    Arguments:
        observable_df: PEtab observable table
        measurement_df: PEtab measurement table
        parameter_df: PEtab measurement table

    Returns:

    """

    observable_parameter_override_id_to_placeholder_id, \
    noise_parameter_override_id_to_placeholder_id = \
        get_parameter_override_id_to_placeholder_id(
            observable_df=observable_df,
            measurement_df=measurement_df)

    # parameters selected for hierarchical optimization
    hierarchical_candidates = parameter_df.index[
        (parameter_df.estimate == 1)
        & (parameter_df.hierarchicalOptimization == 1)]

    offset_candidates = set()
    scaling_candidates = set()
    sigma_candidates = set()

    for optimization_parameter_id in hierarchical_candidates:
        # check which model parameter this one overrides

        if optimization_parameter_id \
                in observable_parameter_override_id_to_placeholder_id:
            placeholder_ids = \
                observable_parameter_override_id_to_placeholder_id[
                    optimization_parameter_id]

            # check in which observables this parameter occurs
            for placeholder_id in placeholder_ids:
                observable_id = '_'.join(placeholder_id.split('_')[1:])
                observable_formula = observable_df.loc[observable_id,
                                                       ptc.OBSERVABLE_FORMULA]

                if parameter_is_offset_parameter(
                        placeholder_id, observable_formula):
                    offset_candidates.add(optimization_parameter_id)
                elif parameter_is_scaling_parameter(
                        placeholder_id, observable_formula):
                    scaling_candidates.add(optimization_parameter_id)
                else:
                    raise RuntimeError(
                        f'Parameter {optimization_parameter_id} selected '
                        'for hierarchical optimization but is neither '
                        'offset, proportionality or sigma parameter in '
                        f'{observable_id}: {observable_formula}.'
                        'Dunno what to do.')
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

    # TODO Can't use hierarchical optimization with non-normal or
    #  transformation yet
    if (offset_candidates or scaling_candidates or sigma_candidates):
        not_normal = ptc.NOISE_DISTRIBUTION in observable_df \
                     and not np.all(x == ptc.NORMAL
                                    or (isinstance(x, Number) and np.isnan(x))
                                    for x
                                    in observable_df[ptc.NOISE_DISTRIBUTION])
        not_lin = ptc.OBSERVABLE_TRANSFORMATION in observable_df \
                     and not np.all(x == ptc.LIN
                                    or (isinstance(x, Number) and np.isnan(x))
                                    for x in observable_df[
                                        ptc.OBSERVABLE_TRANSFORMATION])
        if not_normal or not_lin:
            raise ValueError("Can't use hierarchical optimization with "
                             "non-normal noise or observable transformation "
                             "yet.")

    return (list(offset_candidates),
            list(scaling_candidates),
            list(sigma_candidates))


def get_analytical_parameter_table(
        hierarchical_candidate_ids: list,
        parameter_type: str,
        condition_id_to_index: Dict[str, int],
        measurement_df: pd.DataFrame,
        observable_ids,
        condition_map,
        no_preeq_condition_idx: int
) -> List[Tuple[int, int, int]]:
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

    # need list, not ndarray
    condition_map_list = [list(x) for x in condition_map]

    if parameter_type == 'observable':
        def _get_overrides():
            return split_parameter_replacement_list(row.observableParameters)
    elif parameter_type == 'noise':
        def _get_overrides():
            return split_parameter_replacement_list(row.noiseParameters)
    else:
        raise ValueError("parameter_type must be 'noise' or "
                         f"'observable', but got {parameter_type}")

    use = []
    for _, row in measurement_df.iterrows():
        overrides = _get_overrides()

        sim_cond_idx = \
            condition_id_to_index[row.simulationConditionId]
        preeq_cond_idx = no_preeq_condition_idx
        if not (isinstance(row.preequilibrationConditionId, Number)
                and isnan(row.preequilibrationConditionId)):
            preeq_cond_idx = condition_id_to_index[
                row.preequilibrationConditionId]

        for s in overrides:
            # print(s, parametersForHierarchical)
            try:
                candidate_idx = hierarchical_candidate_ids.index(s)
            except ValueError:
                continue  # current parameter not in list

            condition_idx = condition_map_list.index(
                [preeq_cond_idx, sim_cond_idx])
            observable_idx = observable_ids.index(row.observableId)
            tup = (candidate_idx, condition_idx, observable_idx)

            # Don't add a new line for each timepoint
            # We don't allow separate parameters for individual time-points
            # (Can be implemented via different observables)
            if tup not in use:
                use.append(tup)

    if not len(use):
        raise AssertionError("Candidates were: "
                             f"{hierarchical_candidate_ids} but nothing "
                             "usable found")

    return use
