#!/usr/bin/env python3

import sys
from libsbml import *

def create_unit_definitions(model):
    unitDefinition = model.createUnitDefinition()
    unitDefinition.setId('litre_per_second')
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_LITRE)
    unit.setExponent(1)
    unit.setScale(0)
    unit.setMultiplier(1)
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_SECOND)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)

    unitDefinition = model.createUnitDefinition()
    unitDefinition.setId('mole_per_litre')
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_MOLE)
    unit.setExponent(1)
    unit.setScale(0)
    unit.setMultiplier(1)
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_LITRE)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)


    unitDefinition = model.createUnitDefinition()
    unitDefinition.setId('mole_per_second')
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_MOLE)
    unit.setExponent(1)
    unit.setScale(0)
    unit.setMultiplier(1)
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_SECOND)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)

    unitDefinition = model.createUnitDefinition()
    unitDefinition.setId('litre2_per_mole_per_second')
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_LITRE)
    unit.setExponent(2)
    unit.setScale(0)
    unit.setMultiplier(1)
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_MOLE)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)
    unit = unitDefinition.createUnit()
    unit.setKind(UNIT_KIND_SECOND)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)

def create_species(model, id, compartment, constant = False, initialAmount = 0.0,
                   substanceUnits = 'mole', boundaryCondition = False, hasOnlySubstanceUnits = False):
    s = model.createSpecies()
    s.setId(id)
    s.setCompartment(compartment)
    s.setConstant(constant)
    s.setInitialAmount(initialAmount)
    s.setSubstanceUnits(substanceUnits)
    s.setBoundaryCondition(boundaryCondition)
    s.setHasOnlySubstanceUnits(hasOnlySubstanceUnits)
    return s

def create_parameter(model, id, constant, value, units):
    k = model.createParameter()
    k.setId(id)
    k.setName(id)
    k.setConstant(constant)
    k.setValue(value)
    k.setUnits(units)
    return k

def create_reaction(model, id, reactants, products, formula, reversible = False, fast = False):
    r = model.createReaction()
    r.setId(id)
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

    math_ast = parseL3Formula(formula)
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

def add_sigma(model, observable):
    # corresponding sigma
    p = create_parameter(model, f'sigma_{observable}', False, 1.0,
                         'mole_per_litre')
    p = create_parameter(model, f'noiseParameter1_{observable}', True, 0.2,
                         'mole_per_litre')
    rule = create_assigment_rule(model, f'sigma_{observable}',
                                 f'noiseParameter1_{observable}')


def create_model():
    document = SBMLDocument(3, 1)
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

    # write observables
    for i in range(1, 4):
        # all states fully observed
        observable = f'observable_x{i}'
        p = create_parameter(model, observable, False, 1.0, 'mole_per_litre')
        rule = create_assigment_rule(model, observable, 'x%d' % i)
        add_sigma(model, f'x{i}')

    # Scaled x1
    p = create_parameter(model, 'observableParameter1_x1_scaled', True, 2.0, 'dimensionless')
    p = create_parameter(model, 'observable_x1_scaled', False, 1.0, 'mole_per_litre')
    rule = create_assigment_rule(model, 'observable_x1_scaled', 'observableParameter1_x1_scaled * x1')
    add_sigma(model, 'x1_scaled')

    # x2 with offset
    p = create_parameter(model, 'observableParameter1_x2_offsetted', True, 3.0, 'mole_per_litre')
    p = create_parameter(model, 'observable_x2_offsetted', False, 1.0, 'mole_per_litre')
    rule = create_assigment_rule(model, 'observable_x2_offsetted', 'observableParameter1_x2_offsetted + x2')
    add_sigma(model, 'x2_offsetted')

    # Fully observed with sigma parameter
    # observable
    p = create_parameter(model, 'observable_x1withsigma', False, 1.0,
                         'mole_per_litre')
    rule = create_assigment_rule(model, 'observable_x1withsigma', 'x1')
    add_sigma(model, 'x1withsigma')

    return writeSBMLToString(document)


if __name__ == '__main__':
    print(create_model())
