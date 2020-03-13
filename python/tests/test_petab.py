from petab.C import *
import parpe.petab
import pandas as pd


def test_get_parameter_override_id_to_placeholder_id():
    """Test get_parameter_override_id_to_placeholder_id"""

    observable_df = pd.DataFrame(data={
        OBSERVABLE_ID: ['obs1', 'obs2'],
        OBSERVABLE_FORMULA: ['1.0',
                             'observableParameter1_obs2 + '
                             'observableParameter2_obs2'],
        NOISE_FORMULA: ['noiseParameter1_obs1', '1.0']
    })
    observable_df.set_index(OBSERVABLE_ID, inplace=True)

    measurement_df = pd.DataFrame(data={
        OBSERVABLE_ID: ['obs1', 'obs1', 'obs2', 'obs2'],
        OBSERVABLE_PARAMETERS: ['', '', '1.0;1.0', 'override1;override2'],
        NOISE_PARAMETERS: ['1.0', 'override3', '', '']
    })

    expected = ({'override1': ['observableParameter1_obs2'],
                 'override2': ['observableParameter2_obs2']},
                {'override3': ['noiseParameter1_obs1']})

    actual = parpe.petab.get_parameter_override_id_to_placeholder_id(
        observable_df=observable_df, measurement_df=measurement_df)

    assert expected == actual


def test_get_observable_placeholders_by_observable():
    """Test get_observable_placeholders_by_observable"""
    observable_df = pd.DataFrame(data={
        OBSERVABLE_ID: ['obs1', 'obs2'],
        OBSERVABLE_FORMULA: ['1.0',
                             'observableParameter1_obs2 + '
                             'observableParameter2_obs2'],
        NOISE_FORMULA: ['noiseParameter1_obs1', '1.0']
    })
    observable_df.set_index(OBSERVABLE_ID, inplace=True)

    expected = {'obs1': [],
                'obs2': ['observableParameter1_obs2',
                         'observableParameter2_obs2']}

    actual = parpe.petab.get_observable_placeholders_by_observable(
        observable_df)

    assert actual == expected

    expected = {'obs1': ['noiseParameter1_obs1'],
                'obs2': []}

    actual = parpe.petab.get_noise_placeholders_by_observable(
        observable_df)

    assert actual == expected
