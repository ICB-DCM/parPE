"""Test related to ``petab.hierarchical_optimization``"""

from parpe.hierarchical_optimization import (parameter_is_scaling_parameter,
                                             parameter_is_offset_parameter)


def test_parameter_is_offset_parameter():
    assert parameter_is_offset_parameter('a', 'a + b') is True
    assert parameter_is_offset_parameter('b', 'a + b') is True
    assert parameter_is_offset_parameter('b', 'a - b') is False
    assert parameter_is_offset_parameter('b', 'sqrt(b)') is False
    assert parameter_is_offset_parameter('b', 'a * b') is False


def test_parameter_is_scaling_parameter():
    assert parameter_is_scaling_parameter('a', 'a + b') is False
    assert parameter_is_scaling_parameter('a', 'a * b') is True
    assert parameter_is_scaling_parameter('a', 'a * b + 1') is False
    assert parameter_is_scaling_parameter('a', 'a * a') is False
