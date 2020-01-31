"""Test related to ``petab.hierarchical_optimization``"""

import petab


def test_parameter_is_offset_parameter():
    assert petab.parameter_is_offset_parameter('a', 'a + b') is True
    assert petab.parameter_is_offset_parameter('b', 'a + b') is True
    assert petab.parameter_is_offset_parameter('b', 'a - b') is False
    assert petab.parameter_is_offset_parameter('b', 'sqrt(b)') is False
    assert petab.parameter_is_offset_parameter('b', 'a * b') is False


def test_parameter_is_scaling_parameter():
    assert petab.parameter_is_scaling_parameter('a', 'a + b') is False
    assert petab.parameter_is_scaling_parameter('a', 'a * b') is True
    assert petab.parameter_is_scaling_parameter('a', 'a * b + 1') is False
    assert petab.parameter_is_scaling_parameter('a', 'a * a') is False
