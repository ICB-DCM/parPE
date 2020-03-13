#!/usr/bin/env python3
"""Create SBML model for 'steadystate' example"""

import libsbml


def add_unit(unit_definition: libsbml.UnitDefinition,
             kind, exponent = 1, scale = 0, multiplier = 1) -> libsbml.Unit:
    """Add unit to SBML unit definition"""
    unit = unit_definition.createUnit()
    unit.setKind(kind)
    unit.setExponent(exponent)
    unit.setScale(scale)
    unit.setMultiplier(multiplier)
    return unit


def create_unit_definitions(model: libsbml.Model) -> None:
    """Create SBML unit definitions"""

    unit_definition = model.createUnitDefinition()
    unit_definition.setId('litre_per_second')
    add_unit(unit_definition, kind=libsbml.UNIT_KIND_LITRE)
    add_unit(unit_definition, kind=libsbml.UNIT_KIND_SECOND, exponent=-1)

    unit_definition = model.createUnitDefinition()
    unit_definition.setId('mole_per_litre')
    add_unit(unit_definition, kind=libsbml.UNIT_KIND_MOLE)
    add_unit(unit_definition, kind=libsbml.UNIT_KIND_LITRE, exponent=-1)

    unit_definition = model.createUnitDefinition()
    unit_definition.setId('mole_per_second')
    add_unit(unit_definition, kind=libsbml.UNIT_KIND_MOLE)
    add_unit(unit_definition, kind=libsbml.UNIT_KIND_SECOND, exponent=-1)

    unit_definition = model.createUnitDefinition()
    unit_definition.setId('litre2_per_mole_per_second')
    add_unit(unit_definition, kind=libsbml.UNIT_KIND_LITRE, exponent=2)
    add_unit(unit_definition, kind=libsbml.UNIT_KIND_MOLE, exponent=-1)
    add_unit(unit_definition, kind=libsbml.UNIT_KIND_SECOND, exponent=-1)


def create_species(
        model: libsbml.Model, sbml_id:str, compartment,
        constant = False, initial_amount = 0.0,
        substance_units ='mole', boundary_condition = False,
        has_only_substance_units = False):
    """Create SBML species"""
    s = model.createSpecies()
    s.setId(sbml_id)
    s.setCompartment(compartment)
    s.setConstant(constant)
    s.setInitialAmount(initial_amount)
    s.setSubstanceUnits(substance_units)
    s.setBoundaryCondition(boundary_condition)
    s.setHasOnlySubstanceUnits(has_only_substance_units)
    return s


def create_parameter(model, parameter_id, constant, value, units):
    """Add parameter to SBML model"""
    k = model.createParameter()
    k.setId(parameter_id)
    k.setName(parameter_id)
    k.setConstant(constant)
    k.setValue(value)
    k.setUnits(units)
    return k


def create_reaction(model, reaction_id, reactants, products, formula,
                    reversible = False, fast = False):
    """Add reaction to SBML model"""
    r = model.createReaction()
    r.setId(reaction_id)
    r.setReversible(reversible)
    r.setFast(fast)

    for (coeff, name) in reactants:
        species_ref = r.createReactant()
        species_ref.setSpecies(name)
        species_ref.setConstant(True) # TODO ?
        species_ref.setStoichiometry(coeff)

    for (coeff, name) in products:
        species_ref = r.createProduct()
        species_ref.setSpecies(name)
        species_ref.setConstant(True) # TODO ?
        species_ref.setStoichiometry(coeff)

    math_ast = libsbml.parseL3Formula(formula)
    kinetic_law = r.createKineticLaw()
    kinetic_law.setMath(math_ast)
    return r


def create_assigment_rule(model, name, formula):
    rule = model.createAssignmentRule()
    rule.setId(name)
    rule.setName(name)
    rule.setVariable(name)
    rule.setFormula(formula)
    return rule


def create_model():
    document = libsbml.SBMLDocument(3, 1)
    model = document.createModel()
    model.setTimeUnits("second")
    model.setExtentUnits("mole")
    model.setSubstanceUnits('mole')

    create_unit_definitions(model)

    c1 = model.createCompartment()
    c1.setId('c1')
    c1.setConstant(True)
    c1.setSize(1)
    c1.setSpatialDimensions(3)
    c1.setUnits('litre')

    s1 = create_species(model, 'x1', 'c1', False, 0.1)
    s2 = create_species(model, 'x2', 'c1', False, 0.4)
    s3 = create_species(model, 'x3', 'c1', False, 0.7)

    #TODO: initial amounts should be parameters
    p1 = create_parameter(model, 'p1', True, 1.0, 'litre2_per_mole_per_second')
    p2 = create_parameter(model, 'p2', True, 0.5, 'litre2_per_mole_per_second')
    p3 = create_parameter(model, 'p3', True, 0.4, 'litre_per_second')
    p4 = create_parameter(model, 'p4', True, 2.0, 'litre_per_second')
    p5 = create_parameter(model, 'p5', True, 0.1, 'mole_per_second')
    k0 = create_parameter(model, 'k0', True, 1.0, 'litre_per_second')


    create_reaction(model, 'r1', [(2, 'x1')], [(1, 'x2')], 'p1 * x1^2')
    create_reaction(model, 'r2', [(1, 'x1'), (1, 'x2')], [(1, 'x3')], 'p2 * x1 * x2')
    create_reaction(model, 'r3', [(1, 'x2')], [(2, 'x1')], 'p3 * x2')
    create_reaction(model, 'r4', [(1, 'x3')], [(1, 'x1'), (1, 'x2')], 'p4 * x3')
    create_reaction(model, 'r5', [(1, 'x3')], [], 'k0 * x3')
    create_reaction(model, 'r6', [], [(1, 'x1')], 'p5')

    #    S -2 -1  2  1  0 1  v1: p1 x1^2
    #       1 -1 -1  1  0 0  v2: p2 x1 x2
    #       0  1  0 -1 -1 0  v3: p3 x2
    #                        v4: p4 x3
    #                        v5: k0 x3
    #                        v6: p5
    # R1: 2X1         ->       X2
    # R2:  X1 + X2    ->           X3
    # R3:       X2    -> 2X1
    # R4:          X3 ->  X1 + X2
    # R5:          X3 ->
    # R6:             ->  X1

    return libsbml.writeSBMLToString(document)


if __name__ == '__main__':
    print(create_model())
