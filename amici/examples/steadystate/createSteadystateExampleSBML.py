#!/usr/bin/env python3
 
import sys
from libsbml import *

def create_unit_definitions(model):
    per_second = model.createUnitDefinition()
    per_second.setId('per_second')
    unit = per_second.createUnit()
    unit.setKind(UNIT_KIND_SECOND)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)
    
    mole_per_second = model.createUnitDefinition()
    mole_per_second.setId('mole_per_second')
    unit = mole_per_second.createUnit()
    unit.setKind(UNIT_KIND_MOLE)
    unit.setExponent(1)
    unit.setScale(0)
    unit.setMultiplier(1)
    unit = mole_per_second.createUnit()
    unit.setKind(UNIT_KIND_SECOND)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)
    
    per_mole_per_second = model.createUnitDefinition()
    per_mole_per_second.setId('per_mole_per_second')
    unit = per_mole_per_second.createUnit()
    unit.setKind(UNIT_KIND_MOLE)
    unit.setExponent(-1)
    unit.setScale(0)
    unit.setMultiplier(1)
    unit = per_mole_per_second.createUnit()
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
    rule.setVariable(name)
    rule.setFormula(formula)
    return rule

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
    p1 = create_parameter(model, 'p1', True, 1.0, 'per_mole_per_second')
    p2 = create_parameter(model, 'p2', True, 0.5, 'per_second')
    p3 = create_parameter(model, 'p3', True, 0.4, 'mole_per_second')
    p4 = create_parameter(model, 'p4', True, 2.0, 'per_second')
    p5 = create_parameter(model, 'p5', True, 0.1, 'per_second')
    k4 = create_parameter(model, 'k4', True, 1.0, 'mole_per_second')

    
    create_reaction(model, 'r1', [(2, 'x1')], [(1, 'x2')], 'p1 * x1^2')
    create_reaction(model, 'r2', [(1, 'x1'), (1, 'x2')], [(1, 'x3')], 'p2 * x1 * x2')
    create_reaction(model, 'r3', [(1, 'x2')], [(2, 'x1')], 'p3 * x2')
    create_reaction(model, 'r4', [(1, 'x3')], [(1, 'x1'), (1, 'x2')], 'p4 * x3')
    create_reaction(model, 'r5', [(1, 'x3')], [], 'k4 * x3')
    create_reaction(model, 'r6', [], [(1, 'x1')], 'p5')
    
    #    S -2 -1  2  1  0 1  v1: p1 x1^2 
    #       1 -1 -1  1  0 0  v2: p2 x1
    #       0  1  0 -1 -1 0  v3: p3
    #                        v4: p4 x3
    #                        v5: k4 * x3
    #                        v6: p5
    # R1: 2X1         ->       X2
    # R2:  X1 + X2    ->           X3
    # R3:       X2    -> 2X1 
    # R4:          X3 ->  X1 + X2 
    # R5:          X3 ->
    # R6:             ->  X1

    # write observables
    for i in range(1, 4):
        observable = 'observable_x%d' %i
        p = create_parameter(model, observable, True, 1.0, '')
        rule = create_assigment_rule(model, observable, 'x%d' % i)
    
    p = create_parameter(model, 'scaling_x1', True, 1.0, '')
    p = create_parameter(model, 'observable_x1_scaled', True, 1.0, '')
    rule = create_assigment_rule(model, 'observable_x1_scaled', 'scaling_x1 * x1')
        
    p = create_parameter(model, 'offset_x2', True, 1.0, '')
    p = create_parameter(model, 'observable_x2_offsetted', True, 1.0, '')
    rule = create_assigment_rule(model, 'observable_x2_offsetted', 'offset_x2 + x2')
    
    p = create_parameter(model, 'observable_x1withsigma_sigma', True, 1.0, '')
    p = create_parameter(model, 'observable_x1withsigma', True, 1.0, '')
    rule = create_assigment_rule(model, 'observable_x1withsigma', 'x1')
  
    return writeSBMLToString(document)

if __name__ == '__main__':
    print(create_model())
 