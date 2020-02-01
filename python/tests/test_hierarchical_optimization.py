"""Test related to ``petab.hierarchical_optimization``"""
from numpy import nan
import pandas as pd
from petab import C as ptc
import pytest
from parpe.hierarchical_optimization import (parameter_is_scaling_parameter,
                                             parameter_is_offset_parameter,
                                             get_analytical_parameter_table)


def test_parameter_is_offset_parameter():
    """Test parameter_is_offset_parameter"""

    assert parameter_is_offset_parameter('a', 'a + b') is True
    assert parameter_is_offset_parameter('b', 'a + b') is True
    assert parameter_is_offset_parameter('b', 'a - b') is False
    assert parameter_is_offset_parameter('b', 'sqrt(b)') is False
    assert parameter_is_offset_parameter('b', 'a * b') is False


def test_parameter_is_scaling_parameter():
    """Test parameter_is_scaling_parameter"""

    assert parameter_is_scaling_parameter('a', 'a + b') is False
    assert parameter_is_scaling_parameter('a', 'a * b') is True
    assert parameter_is_scaling_parameter('a', 'a * b + 1') is False
    assert parameter_is_scaling_parameter('a', 'a * a') is False


def test_get_analytical_parameter_table():
    """Test get_analytical_parameter_table"""

    condition_id_to_index = {'c0': 0,
                             'c1': 1,
                             'c2': 2,
                             'c3': 3}
    measurement_df = pd.DataFrame(data={
        ptc.OBSERVABLE_ID: ['x2', 'x2'],
        ptc.SIMULATION_CONDITION_ID: ['c1', 'c2'],
        ptc.PREEQUILIBRATION_CONDITION_ID: [nan, nan],
        ptc.OBSERVABLE_PARAMETERS: ['scaling_x2', 'scaling_x2;1.0'],
        ptc.NOISE_PARAMETERS: ['noise', ''],

    })
    observable_ids = ('x1', 'x2', 'x3')
    condition_map = [[-1, 0],
                     [-1, 1],
                     [-1, 2],
                     [-1, 3]]
    expected = [
        # hierarchical_par_idx, condition_idx, observable_idx
        (0, 1, 1),
        (0, 2, 1),
    ]

    actual = get_analytical_parameter_table(
        hierarchical_candidate_ids=['scaling_x2'],
        parameter_type="observable",
        condition_id_to_index=condition_id_to_index,
        measurement_df=measurement_df,
        observable_ids=observable_ids,
        condition_map=condition_map,
        no_preeq_condition_idx=-1
    )

    assert actual == expected

    # For sigmas
    expected = [
        # hierarchical_par_idx, condition_idx, observable_idx
        (0, 1, 1),
    ]
    actual = get_analytical_parameter_table(
        hierarchical_candidate_ids=['noise'],
        parameter_type="noise",
        condition_id_to_index=condition_id_to_index,
        measurement_df=measurement_df,
        observable_ids=observable_ids,
        condition_map=condition_map,
        no_preeq_condition_idx=-1
    )

    assert actual == expected

    # Test unused parameter
    with pytest.raises(AssertionError):
        get_analytical_parameter_table(
            hierarchical_candidate_ids=['unused_par'],
            parameter_type="observable",
            condition_id_to_index=condition_id_to_index,
            measurement_df=measurement_df,
            observable_ids=observable_ids,
            condition_map=condition_map,
            no_preeq_condition_idx=-1
        )
