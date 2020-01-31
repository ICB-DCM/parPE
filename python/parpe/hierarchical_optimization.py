"""Functions related to hierarchical optimization

https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz581/5538985
"""

import sympy as sp


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

